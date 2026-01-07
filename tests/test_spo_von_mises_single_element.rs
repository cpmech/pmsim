use gemlab::mesh::Samples;
use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_sparse::Genie;
use serial_test::serial;

// IMPORTANT:
// Since MUMPS is not thread-safe, we need to use serial_test::serial

// von Mises plasticity with a single-element
//
// This test runs a plane-strain compression of a single element represented
// by the von Mises model. The results are compared with the code HYPLAS
// discussed in Ref #1.
//
// TEST GOAL
//
// Verifies the plane-strain implementation of the von Mises model.
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
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME: &str = "spo_von_mises_single_element";

// constants
const YOUNG: f64 = 1500.0;
const POISSON: f64 = 0.25;
const Z_INI: f64 = 9.0;
const NU: f64 = POISSON;
const NU2: f64 = POISSON * POISSON;
const NGAUSS: usize = 1;
const NSTAGE: usize = 5;

#[test]
#[serial]
fn test_spo_von_mises_single_element() -> Result<(), StrError> {
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
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // stage-wise vertical displacement increment
    let delta_y = -Z_INI * (1.0 - NU2) / (YOUNG * f64::sqrt(1.0 - NU + NU2));

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&top, Dof::Uy, |t| delta_y * t);

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_lagrange_mult_method(true)
        .set_steady(NSTAGE)
        .set_max_iterations(20);

    // solve and check with UMFPACK
    solve_and_check(&mesh, &schema, &essential, &natural, &config)?;

    // solve and check with MUMPS
    config.set_lin_sol_genie(Genie::Mumps);
    solve_and_check(&mesh, &schema, &essential, &natural, &config)?;
    Ok(())
}

fn solve_and_check(
    mesh: &Mesh,
    schema: &Schema,
    essential: &BcEssential,
    natural: &BcNatural,
    config: &Config,
) -> Result<(), StrError> {
    // solution
    SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?;

    // compare the results with Ref #1
    let tol_displacement = 1e-13;
    let tol_stress = 1e-10;
    let all_good = compare_results(
        &mesh,
        &schema,
        &config,
        "/tmp/pmsim/",
        NAME,
        ReferenceDataType::SPO,
        &format!("data/spo/{}.json", NAME),
        tol_displacement,
        tol_stress,
        0,
    )?;
    assert!(all_good);
    Ok(())
}
