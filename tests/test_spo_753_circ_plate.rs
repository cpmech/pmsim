use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::compare_results;
use russell_lab::*;

const NAME: &str = "test_spo_753_circ_plate";
const SAVE_FIGURE: bool = false;

const P_ARRAY: [f64; 13] = [
    0.0, 100.0, 200.0, 220.0, 230.0, 240.0, 250.0, 255.0, 257.0, 259.0, 259.5, 259.75, 259.77,
];
const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.24; // Poisson's coefficient
const Z_INI: f64 = 16000.0; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_753_circ_plate() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/meshes/{}.msh", NAME))?;
    if SAVE_FIGURE {
        mesh.check_all()?;
        let mut fig = Figure::new();
        fig.size(1000.0, 1000.0)
            .draw(&mesh, &format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    // features
    let features = Features::new(&mesh, false);
    let top = features.search_edges(At::Y(1.0), any_x)?;
    let right_corner = features.search_point_ids(At::XY(10.0, 0.0), any_x)?;
    assert_eq!(right_corner.len(), 1);

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.point(right_corner[0], Dof::Uy, 0.0);
    // println!("{}", essential);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges_fn(&top, Nbc::Qn, |t| -P_ARRAY[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        .set_incremental(P_ARRAY.len())
        .set_lagrange_mult_method(false)
        .set_tol_rr(1e-6)
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, NAME, None)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    return Ok(());

    // verify the results
    let tol_displacement = 1e-1;
    let tol_stress = 1e-1;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        "spo_753_circ_plate.json",
        tol_displacement,
        tol_stress,
        2,
    )?;
    assert!(all_good);
    Ok(())
}
