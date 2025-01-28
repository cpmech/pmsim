use gemlab::prelude::*;
use pmsim::base::{Config, Dof, Elem, Essential, Natural, Nbc, ParamSolid, StressStrain};
use pmsim::fem::{BcDistributedArray, BcPrescribedArray, Elements, FemMesh, FemState, LinearSystem};
use russell_lab::*;
use russell_sparse::prelude::*;

const OUT_DIR: &str = "/tmp/pmsim/matrix-market";
const SAVE_FIGURE: bool = false;

const R1: f64 = 3.0; // inner radius
const R2: f64 = 6.0; // outer radius
const P1: f64 = 200.0; // inner pressure (magnitude)
const P2: f64 = 100.0; // outer pressure (magnitude)
const YOUNG: f64 = 1000.0; // Young's modulus
const POISSON: f64 = 0.25; // Poisson's coefficient
const NA: usize = 94; // number of alpha divisions

fn generate_matrix(name: &str, nr: usize) -> Result<SparseMatrix, StrError> {
    // generate mesh
    let mesh = Structured::quarter_ring_2d(R1, R2, nr, NA, GeoKind::Qua4, true).unwrap();

    // draw mesh
    if SAVE_FIGURE {
        let mut fig = Figure::new();
        fig.figure_size = Some((800.0, 800.0));
        mesh.draw(Some(fig), &format!("{}/{}.svg", OUT_DIR, name), |_, _| {})?;
    }

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, R1), any_x)?;
    let outer_circle = features.search_edges(At::Circle(0.0, 0.0, R2), any_x)?;

    // material parameters
    let param1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
    };
    let fem = FemMesh::new(&mesh, [(1, Elem::Solid(param1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // prescribed values
    let prescribed_values = BcPrescribedArray::new(&fem, &essential)?;

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .edges(&inner_circle, Nbc::Qn, -P1)
        .edges(&outer_circle, Nbc::Qn, -P2);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lin_sol_genie(Genie::Umfpack)
        .access_lin_sol_params()
        .umfpack_enforce_unsymmetric_strategy = true;

    // elements
    let mut elements = Elements::new(&fem, &config)?;

    // boundaries
    let mut boundaries = BcDistributedArray::new(&fem, &config, &natural)?;

    // FEM state
    let state = FemState::new(&fem, &config)?;

    // compute jacobians
    elements.calc_jacobians(&state)?;
    boundaries.calc_jacobians(&state)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&fem, &config, &prescribed_values, &elements, &boundaries)?;

    // assemble jacobian matrix
    let kk = lin_sys.jacobian.get_coo_mut()?;
    elements.assemble_jacobians(kk, &prescribed_values.flags)?;
    boundaries.assemble_jacobians(kk, &prescribed_values.flags)?;

    // augment global Jacobian matrix
    for eq in &prescribed_values.equations {
        kk.put(*eq, *eq, 1.0).unwrap();
    }

    // done
    Ok(lin_sys.jacobian)
}

fn run(name: &str, mat: &mut SparseMatrix, enforce_unsym_strategy: bool) -> Result<(), StrError> {
    // allocate stats structure
    let mut stats = StatsLinSol::new();

    // save information about the matrix
    let (nrow, ncol, nnz, sym) = mat.get_info();
    stats.matrix.name = name.to_string();
    stats.matrix.nrow = nrow;
    stats.matrix.ncol = ncol;
    stats.matrix.nnz = nnz;
    stats.matrix.symmetric = format!("{:?}", sym);

    // lin solver params
    let mut params = LinSolParams::new();
    params.umfpack_enforce_unsymmetric_strategy = enforce_unsym_strategy;

    // allocate and configure the solver
    let mut solver = LinSolver::new(Genie::Umfpack)?;

    // call factorize
    solver.actual.factorize(mat, Some(params))?;

    // allocate vectors
    let mut x = Vector::new(nrow);
    let rhs = Vector::filled(nrow, 1.0);

    // solve linear system
    solver.actual.solve(&mut x, mat, &rhs, false)?;

    // verify the solution
    stats.verify = VerifyLinSys::from(mat, &x, &rhs)?;

    // update and print stats
    solver.actual.update_stats(&mut stats);
    println!("{}", stats.get_json());
    let path = if enforce_unsym_strategy {
        format!("/tmp/pmsim/{}-enforce-unsym.json", name)
    } else {
        format!("/tmp/pmsim/{}-auto.json", name)
    };
    stats.write_json(&path).unwrap();

    // check
    if enforce_unsym_strategy {
        println!("\ndiff = {:e}\n", stats.verify.max_abs_diff);
        if name == "pres-cylin-bad" {
            assert!(stats.verify.max_abs_diff > 4100.0);
        } else {
            assert!(stats.verify.max_abs_diff < 1e-4);
        }
    } else {
        if name == "pres-cylin-bad" {
            assert!(stats.verify.max_abs_diff < 1e-11);
        } else {
            assert!(stats.verify.max_abs_diff < 1e-11);
        }
    }
    Ok(())
}

fn main() -> Result<(), StrError> {
    let nrs = [
        ("pres-cylin-good", 54), //
        ("pres-cylin-bad", 55),  //
    ];
    for (name, nr) in nrs {
        // generate matrix
        let mut mat = generate_matrix(name, nr)?;

        // solve linear system
        run(name, &mut mat, false)?;
        run(name, &mut mat, true)?;

        // write smat and mtx files
        let csc = mat.get_csc()?;
        csc.write_matrix_market(&format!("{}/{}.mtx", OUT_DIR, name), false)?;
        csc.write_matrix_market(&format!("{}/{}.smat", OUT_DIR, name), true)?;
    }
    Ok(())
}
