use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::material::{Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::StrError;

// von Mises plasticity with a single-element
//
// This test runs a plane-strain compression of a single element represented
// by the von Mises model.
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model.
// Also verifies the output of results.
//
// MESH
//
// Unit square
//
// displacement    displacement
//         ↓         ↓
//  roller 3---------2
//         |         |   E = 1500  z0 = 9.0
//         |         |   ν = 0.25  H = 800
//         |         |
//         0---------1
//      fixed       roller
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

const NAME: &str = "test_von_mises_single_element_2d";
const SAVE_FIGURE: bool = true;

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 1;
const NSTAGE: usize = 5;

#[test]
fn test_von_mises_single_element_2d() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(1.0), any_x)?;

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
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

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
        .set_out_local_state(0)
        .set_lagrange_mult_method(true)
        .set_steady(NSTAGE)
        .set_substepping(true)
        .set_max_iterations(20)
        .update_model_settings(1)
        .set_save_strain(true);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &base, &config)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // analysis
    let ss = results.get_local_state(0).unwrap();
    for s in ss {
        println!("{:?}", s.elastic);
    }

    // figure
    if SAVE_FIGURE {
        let data = PlotterData::from_states(ss);
        let mut plotter = Plotter::new();
        plotter.add_3x2(&data, false, |_, _, _| {})?;
        plotter.save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}
