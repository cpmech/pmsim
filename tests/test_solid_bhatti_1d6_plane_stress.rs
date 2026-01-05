use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{vec_approx_eq, Vector};

// Bhatti's Example 1.6 on page 32
//
// Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
//
// TEST GOAL
//
// This test verifies the equilibrium of a thin bracket modelled by assuming plane-stress
//
// MESH
//
// 2.0  fixed 1'-,_load                connectivity:
//            |     '-,_      load      eid : vertices
// 1.5 - - -  |        ,'3-,__            0 :  0, 2, 3
//            |  1   ,'  |    '-,_        1 :  3, 1, 0
// 1.0 - - -  |    ,'    |  3   ,-'5      2 :  2, 4, 5
//            |  ,'  0   |   ,-'   |      3 :  5, 3, 2
//            |,'        |,-'   2  |
// 0.0  fixed 0----------2---------4   constraints:
//           0.0        2.0       4.0   fixed on x and y
//
// BOUNDARY CONDITIONS
//
// Fully fixed @ points 0 and 1
// Distributed load along edges (1,3) and (3,5) with Qn = -20
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 10,000
// Poisson = 0.2
// Plane-stress with thickness = 0.25

#[test]
fn test_solid_bhatti_1d6_plane_stress() -> Result<(), StrError> {
    println!("\n################################### OLD SOLVER ###################################\n");
    run_test(false, false, false)?;
    println!("\n##################################### NATURAL ####################################\n");
    run_test(true, false, false)?; // Natural continuation
    println!("\n################################ ARCLENGTH FULL ##################################\n");
    run_test(true, true, false)?; // Pseudo-arclength continuation without bordering
    println!("\n############################# ARCLENGTH BORDERING ################################\n");
    run_test(true, true, true)?; // Pseudo-arclength continuation with bordering
    Ok(())
}

fn run_test(new_solver: bool, arclength: bool, bordering: bool) -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_1d6_bracket();
    let features = Features::new(&mesh, false);
    let top = Edges {
        all: vec![features.get_edge(1, 3), features.get_edge(3, 5)],
    };

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 10_000.0,
            poisson: 0.2,
        },
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.points(&[0, 1], Dof::Ux, 0.0).points(&[0, 1], Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = BcNatural::new();
    natural.edges(&top, Nbc::Qn, -20.0);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_lagrange_mult_method(true).set_plane_stress(0.25);

    // nonlinear solver configuration
    config.nl_config().set_verbose(true, true, false);
    if arclength {
        config
            .nl_config()
            .set_method(NlMethod::Arclength)
            .set_bordering(bordering);
        if !bordering {
            config.set_ignore_symmetry(true);
        }
    };

    // solution
    let u = if new_solver {
        let state = solve(&mesh, &schema, &config, &essential, &natural)?;
        Vector::from(&&state.u.as_data()[..12])
    } else {
        let mut state = FemState::new(&mesh, &schema, &essential, &config)?;
        let mut results = FemResults::new(&mesh, &schema, &config)?;
        let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
        solver.solve_sys(&mut state, &mut results)?;
        Vector::from(&&state.u.as_data()[..12])
    };

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
        -1.035527877607004e-02, -2.552969847657423e-02,
         4.727650463081949e-03, -2.473565538172127e-02,
        -1.313941349422282e-02, -5.549310752960183e-02,
         8.389015766816341e-05, -5.556637423271112e-02
    ];
    vec_approx_eq(&u, uu_correct, 1e-15);
    Ok(())
}
