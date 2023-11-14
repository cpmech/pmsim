use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use russell_lab::*;

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
// Flux Qt = 8,000 on left side, edge (0,10,11)
// Convection Cc = (55, 20) on top edges (0,2,1), (2,4,3), and (4,6,5)
// Prescribed temperature T = 110 on the bottom edge (8,10,9)
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// Source = 5e6 over the region
// Constant conductivity kx = ky = 45

fn main() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_6d22_heat();
    let feat = Features::new(&mesh, false); // boundary only
    let bottom = feat.search_edges(At::Y(0.0), any_x)?;
    let edges_flux = feat.search_edges(At::X(0.0), any_x)?;
    let edges_conv = vec![
        feat.search_edges(At::Y(0.03), any_x)?.as_slice(),  // top-horizontal
        feat.search_edges(At::X(0.03), any_x)?.as_slice(),  // middle-vertical
        feat.search_edges(At::Y(0.015), any_x)?.as_slice(), // middle-horizontal
    ]
    .concat();

    // input data
    let (kx, ky) = (45.0, 45.0);
    let source = 5e6;
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::Constant { kx, ky, kz: 0.0 },
        source: Some(source),
    };
    let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&bottom, Ebc::T(|_| 110.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .on(&edges_flux, Nbc::Qt(|_| 8000.0))
        .on(&edges_conv, Nbc::Cv(55.0, |_| 20.0));

    // configuration
    let config = Config::new();

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check U vector
    let tt_bhatti = Vector::from(&[
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
    ]);
    vec_approx_eq(state.uu.as_data(), tt_bhatti.as_data(), 1e-12);
    Ok(())
}
