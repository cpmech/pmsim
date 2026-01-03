use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;

// von Mises plasticity with a four Qua8 elements
//
// This test runs a plane-strain compression of a square represented
// by the von Mises model. The results are compared with the code HYPLAS
// discussed in Ref #1.
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model.
//
// MESH
//
//                 prescribed vertical displacement
//                 ↓       ↓       ↓       ↓       ↓
// 2.0   fix ux > 14------16------13------20------18
//                 |               |               |
//                 |               |               |
// 1.5   fix ux > 17      [2]     15      [3]     19
//                 |               |               |
//                 |               |               |
// 1.0   fix ux >  3-------6-------2------12-------9
//                 |               |               |
//                 |               |               |
// 0.5   fix ux >  7      [0]      5      [1]     11
//                 |               |               |
//                 |               |               |
// 0.0   fix ux >  0-------4-------1------10-------8
//                 ^       ^       ^       ^       ^
//             fix uy  fix uy  fix uy  fix uy  fix uy
//
//                0.0     0.5     1.0     1.5     2.0
//
// xmin = 0.0, xmax = 2.0          E = 1500  z0 = 9.0
// ymin = 0.0, ymax = 2.0          ν = 0.25  H = 800
//
// BOUNDARY CONDITIONS
//
// * Vertically restrain the bottom edge
// * Horizontally restrain the left edge
// * Apply a vertical displacement -δy on the top edge
// * δy is computed such that the first loading will
//   bring the stress point to the yield surface
//
// CONFIGURATION AND PARAMETERS
//
// * Static non-linear plane-strain simulation
// * Young: E = 1500, Poisson: ν = 0.25
// * Hardening: H = 800, Initial yield stress: z0 = 9.0
//
// The results are compared with the code HYPLAS discussed in Ref #1.
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME: &str = "spo_von_mises_2x2_elements";

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 4;
const NSTAGE: usize = 5;

#[test]
fn test_spo_von_mises_2x2_elements() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::block_2d_four_qua8();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(2.0), any_x)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: 800.0,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // stage-wise vertical displacement increment
    let delta_y = -Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&top, Dof::Uy, 1.0, |s, _| delta_y * ((1 + s) as f64));

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_lagrange_mult_method(true)
        .set_steady(NSTAGE)
        .set_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &schema, &config)?;

    // solution
    let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // compare the results with Ref #1
    let tol_displacement = 1e-12;
    let tol_stress = 1e-9;
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        &format!("/tmp/pmsim/{}.json", NAME),
        ReferenceDataType::SPO,
        &format!("data/spo/{}.json", NAME),
        tol_displacement,
        tol_stress,
        0,
    )?;
    assert!(all_good);
    Ok(())
}
