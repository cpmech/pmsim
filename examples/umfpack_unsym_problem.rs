use gemlab::prelude::*;
use pmsim::base::{Config, Ebc, Element, Essential, Natural, Nbc, ParamSolid, ParamStressStrain};
use pmsim::fem::{Boundaries, Data, Elements, LinearSystem, PrescribedValues, State};
use russell_lab::*;
use russell_sparse::prelude::*;

const OUT_DIR: &str = "/tmp/pmsim/matrix-market";
const R1: f64 = 3.0; // inner radius
const R2: f64 = 6.0; // outer radius
const P1: f64 = 200.0; // inner pressure (magnitude)
const P2: f64 = 100.0; // outer pressure (magnitude)
const YOUNG: f64 = 1000.0; // Young's modulus
const POISSON: f64 = 0.25; // Poisson's coefficient
const NA: usize = 94; // number of alpha divisions

fn generate_matrix(nr: usize) -> Result<SparseMatrix, StrError> {
    // generate mesh
    let mesh = Structured::quarter_ring_2d(R1, R2, nr, NA, GeoKind::Qua4).unwrap();

    // features
    let feat = Features::new(&mesh, false);
    let bottom = feat.search_edges(At::Y(0.0), any_x)?;
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let inner_circle = feat.search_edges(At::Circle(0.0, 0.0, R1), any_x)?;
    let outer_circle = feat.search_edges(At::Circle(0.0, 0.0, R2), any_x)?;

    // material parameters
    let param1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(param1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0));

    // prescribed values
    let prescribed_values = PrescribedValues::new(&data, &essential)?;

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .on(&inner_circle, Nbc::Qn(|_| -P1))
        .on(&outer_circle, Nbc::Qn(|_| -P2));

    // configuration
    let mut config = Config::new();
    config.lin_sol_genie = Genie::Umfpack;
    config.lin_sol_params.umfpack_enforce_unsymmetric_strategy = true;

    // elements
    let mut elements = Elements::new(&data, &config)?;

    // boundaries
    let mut boundaries = Boundaries::new(&data, &config, &natural)?;

    // simulation state
    let state = State::new(&data, &config)?;

    // compute jacobians in parallel
    elements.calc_jacobians_parallel(&state)?;
    boundaries.calc_jacobians_parallel(&state)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &config, &prescribed_values, &elements, &boundaries)?;

    // assemble jacobian matrix
    let kk = lin_sys.jacobian.get_coo_mut()?;
    elements.assemble_jacobians(kk, &prescribed_values.flags);
    boundaries.assemble_jacobians(kk, &prescribed_values.flags);

    // augment global Jacobian matrix
    for eq in &prescribed_values.equations {
        kk.put(*eq, *eq, 1.0).unwrap();
    }

    // done
    Ok(lin_sys.jacobian)
}

fn run(name: &str, nr: usize) -> Result<(), StrError> {
    // generate matrix
    let mut mat = generate_matrix(nr)?;

    // allocate stats structure
    let mut stats = StatsLinSol::new();

    // save information about the matrix
    let (nrow, ncol, nnz, symmetry) = mat.get_info();
    stats.matrix.name = name.to_string();
    stats.matrix.nrow = nrow;
    stats.matrix.ncol = ncol;
    stats.matrix.nnz = nnz;
    stats.matrix.symmetry = format!("{:?}", symmetry);

    // lin solver params
    let mut params = LinSolParams::new();
    params.umfpack_enforce_unsymmetric_strategy = true;

    // allocate and configure the solver
    let mut solver = LinSolver::new(Genie::Umfpack)?;

    // call factorize
    solver.actual.factorize(&mut mat, Some(params))?;

    // allocate vectors
    let mut x = Vector::new(nrow);
    let rhs = Vector::filled(nrow, 1.0);

    // solve linear system
    solver.actual.solve(&mut x, &mat, &rhs, false)?;

    // verify the solution
    stats.verify = VerifyLinSys::new(&mat, &x, &rhs)?;

    // update and print stats
    solver.actual.update_stats(&mut stats);
    println!("{}", stats.get_json());

    // write smat and mtx files
    let csr = mat.get_csr_or_from_coo()?;
    csr.write_matrix_market(&format!("{}/{}.mtx", OUT_DIR, name), false)?;

    // check
    if name == "pres-cylin-bad" {
        assert!(stats.verify.max_abs_diff > 1600.0);
    } else {
        assert!(stats.verify.max_abs_diff < 1e-6);
    }
    Ok(())
}

fn main() -> Result<(), StrError> {
    let nrs = [("pres-cylin-good", 54), ("pres-cylin-bad", 55)];
    for (name, nr) in nrs {
        run(name, nr)?;
    }
    Ok(())
}
