use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

// Smith's Example 5.7 (Figure 5.7) on page 178
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a plane-strain simulation with Tri15 elements
//
// MESH
//
//         1.0 kN/m²
//          ↓↓↓↓↓↓
//  0.0  Ux o----o---------------o Ux
//       F  |   /|           _.-'| F
//       I  |  / |       _.-'    | I    15-node
//       X  | /  |   _.-'        | X    triangles
//       E  |/   |.-'            | E
// -2.0  D  o----o---------------o D
//         0.0  1.0             6.0
//             Ux and Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix right edge horizontally
// Fix bottom edge horizontally and vertically
// Concentrated load (Fy) on points 0, 5, 10, 15, 20 equal to
// -0.0778, -0.3556, -0.1333, -0.3556, -0.0778, respectively
//
// NOTE: the distributed load is directly modelled by concentrated forces
// just so we can compare the numeric results with the book results.
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1e5
// Poisson = 0.2
// Plane-strain
//
// NOTE: the Poisson coefficient in the book's figure is different than the
// coefficient in the code. The results given in the book's Fig 5.8 correspond
// to the code's coefficient (Poisson = 0.2)

#[test]
fn test_solid_smith_5d7_tri15_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d7_tri15();

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let right = feat.search_edges(At::X(6.0), any_x)?;
    let bottom = feat.search_edges(At::Y(-2.0), any_x)?;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 1e5,
            poisson: 0.2,
        },
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(|_| 0.0))
        .on(&right, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .at(&[0, 20], Pbc::Fy(|_| -0.0778))
        .at(&[5, 15], Pbc::Fy(|_| -0.3556))
        .at(&[10], Pbc::Fy(|_| -0.1333));

    // configuration
    let mut config = Config::new();
    config.n_integ_point.insert(1, 12);

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00, -1.591418422482600e-05,
         0.000000000000000e+00, -1.158192432406090e-05,
         0.000000000000000e+00, -7.226219651004581e-06,
         0.000000000000000e+00, -3.353677798367540e-06,
         0.000000000000000e+00,  0.000000000000000e+00,
        -9.320588833256607e-07, -1.559179627973348e-05,
         1.493185883616356e-07, -1.127682342509484e-05,
         4.539740902110330e-07, -7.019142981881202e-06,
         3.347179437312136e-07, -3.254591177934953e-06,
         0.000000000000000e+00,  0.000000000000000e+00,
        -1.888048914904935e-06, -1.463414936879944e-05,
         3.230785034937527e-07, -1.040083652084754e-05,
         8.541928153732595e-07, -6.503624127862316e-06,
         6.511070399218456e-07, -3.028030757826080e-06,
         0.000000000000000e+00,  0.000000000000000e+00,
        -2.688953235836881e-06, -1.259279725572649e-05,
         4.988120788393424e-07, -8.864470537186998e-06,
         1.188825956522707e-06, -5.660173178340556e-06,
         9.050868671538049e-07, -2.698165058507540e-06,
         0.000000000000000e+00,  0.000000000000000e+00,
        -2.962033979059817e-06, -8.437865325134800e-06,
         5.008524096747439e-07, -6.728658349863283e-06,
         1.370882667365673e-06, -4.643430624368091e-06,
         1.075230783578103e-06, -2.332827471750362e-06,
         0.000000000000000e+00,  0.000000000000000e+00,
        -1.163757137316748e-06, -4.883879684132889e-07,
         9.372890762240453e-08, -6.201050720099007e-07,
         6.986804775604292e-07, -6.132629410667546e-07,
         6.449990663098843e-07, -4.334700213074987e-07,
         0.000000000000000e+00,  0.000000000000000e+00,
        -1.637918830054012e-07,  2.388485714301134e-07,
        -5.386599282659269e-08,  1.893026527323383e-07,
         1.279441447957562e-07,  1.101243922171629e-07,
         1.655523197871977e-07,  2.727036635857642e-08,
         0.000000000000000e+00,  0.000000000000000e+00,
         1.471290326118162e-07, -9.032439353694393e-08,
         2.512822743294484e-08,  8.102928273256259e-08,
        -3.082728217534403e-09,  7.027436563697669e-08,
        -4.914429394593262e-09,  3.574661176707787e-08,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  4.021675885881946e-08,
         0.000000000000000e+00, -3.259554374873301e-08,
         0.000000000000000e+00, -3.930567433422242e-08,
         0.000000000000000e+00, -1.905778315565390e-08,
         0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-11);
    Ok(())
}
