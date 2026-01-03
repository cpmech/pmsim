use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;

const NAME: &str = "spo_754_footing";
const DRAW_MESH_AND_EXIT: bool = false;
const VERBOSE_LEVEL: usize = 0;

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

// displacement control
const UY: [f64; 14] = [
    -0.01,  // stage =  0
    -0.015, // stage =  1
    -0.02,  // stage =  2
    -0.025, // stage =  3
    -0.035, // stage =  4
    -0.045, // stage =  5
    -0.055, // stage =  6
    -0.065, // stage =  7
    -0.075, // stage =  8
    -0.08,  // stage =  9
    -0.09,  // stage = 10
    -0.11,  // stage = 11
    -0.14,  // stage = 12
    -0.2,   // stage = 13
];

#[test]
fn test_spo_754_footing() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;
    if DRAW_MESH_AND_EXIT {
        mesh.check_all()?;
        let mut draw = Draw::new();
        return draw
            .set_size(800.0, 800.0)
            .zoom_2d(15.0, 69.0, 448.0, 502.0, 0.5, 0.5, 0.5, 0.5)
            .set_range_2d(-10.0, 600.0, -10.0, 600.0)
            .all(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME));
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(500.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let footing = features.search_edges(At::Y(500.0), |x| x[0] <= 50.0)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            // stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&footing, Dof::Uy, 1.0, |stage, _| UY[stage]);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_lagrange_mult_method(true)
        .set_steady(UY.len())
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &schema, &config)?;

    // solution
    let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // verify the results
    let tol_displacement = 1e-10;
    let tol_stress = 5e-5;
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        &format!("/tmp/pmsim/{}.json", NAME),
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);
    Ok(())
}
