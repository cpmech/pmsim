use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::array_approx_eq;

// Bhatti's Example 6.22 on page 449
//
// Bhatti, M.A. (2005) Fundamental Finite Element Analysis and Applications, Wiley, 700p.
//
// TEST GOAL
//
// This test verifies the steady heat equation with prescribed temperature, convection,
// flux, and a volumetric source term. Also, it checks the use of Qua8 elements.
//
// MESH
//
//       0.0    0.015    0.03
// 0.03   0-------1-------2
//        |               |
//        |               3
//        |               |
// 0.015 11            _.'4-------5-------6 0.015
//        |        _.-'                   |
//        |    _.-12                      7 0.0075
//        |_.-'                           |
// 0.0   10---------------9---------------8 0.0
//       0.0             0.03            0.06
//
// BOUNDARY CONDITIONS (see page 445)
//
// Inward flux Qt = -8,000 on left side, edge (0,10,11)
// Convection Cc = (55, 20) on top edges (0,2,1), (2,4,3), and (4,6,5)
// Prescribed temperature T = 110 on the bottom edge (8,10,9)
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// Source = 5e6 over the region
// Constant conductivity kx = ky = 45

#[test]
fn test_heat_bhatti_6d22_convection_sim() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_6d22_heat();
    let features = Features::new(&mesh, false); // boundary only
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let edges_flux = features.search_edges(At::X(0.0), any_x)?;
    let edges_conv_a = features.search_edges(At::Y(0.03), any_x)?; // top-horizontal
    let edges_conv_b = features.search_edges(At::X(0.03), any_x)?; // middle-vertical
    let edges_conv_c = features.search_edges(At::Y(0.015), any_x)?; // middle-horizontal

    // parameters
    let (kx, ky) = (45.0, 45.0);
    let source = 5e6;
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: Some(source),
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;
    let mut config = Config::new(&mesh);
    config.set_lagrange_mult_method(true);

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&bottom, Dof::Phi, 110.0);

    // natural boundary conditions
    let mut natural = BcNatural::new();
    natural
        .edges(&edges_flux, Nbc::Qt, -8000.0) // negative values means inward flux
        .edges(&edges_conv_a, Nbc::Cv(55.0), 20.0)
        .edges(&edges_conv_b, Nbc::Cv(55.0), 20.0)
        .edges(&edges_conv_c, Nbc::Cv(55.0), 20.0);

    // FEM state
    let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

    // FEM results
    let mut results = OutputFiles::new(&mesh, &schema, &config)?;

    // solution
    let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
    solver.solve_sys(&mut state, &mut results)?;

    // check U vector
    let tt_bhatti = &[
        156.440502466202,
        150.75605418729847,
        149.19646294563637,
        144.2245542836661,
        133.8432701060946,
        124.00195294431063,
        121.74635727622194,
        119.14813150652589,
        110.0,
        110.0,
        110.0,
        144.67542222443012,
        129.13200798820264,
    ];
    array_approx_eq(&state.u.as_data()[..13], tt_bhatti, 1e-12);
    Ok(())
}
