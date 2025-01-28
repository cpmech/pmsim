use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

// Smith's Example 5.15 (Figure 5.15) on page 183
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a plane-strain simulation with Qua8 elements
// and reduced integration.
//
// MESH
//
//          1.0 kN/m²
//         ↓↓↓↓↓↓↓↓↓↓↓
//  0.0    0----1----2----3----4
//         |         |         |
//         5         6         7
//         |         |         |
// -3.0 Ux 8----9---10---11---12 Ux
//      F  |         |         | F
//      I 13        14        15 I
//      X  |         |         | X
// -6.0 E 16---17---18---19---20 E
//      D  |         |         | D
//        21        22        23
//         |         |         |
// -9.0   24---25---26---27---28
//        0.0       3.0       6.0
//            Ux and Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix right edge horizontally
// Fix bottom edge horizontally and vertically
// Distributed load Qn = -1 on top edge with x ≤ 3
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1e6
// Poisson = 0.3
// Plane-strain
// NOTE: using reduced integration with 4 points

#[test]
fn test_solid_smith_5d15_qua8_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d15_qua8();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(6.0), any_x)?;
    let bottom = features.search_edges(At::Y(-9.0), any_x)?;
    let top = features.search_edges(At::Y(0.0), |x| x[0] <= 3.0)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 1e6,
            poisson: 0.3,
        },
    };
    let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&top, Nbc::Qn, -1.0);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_ngauss(1, 4);

    // FEM state
    let mut state = FemState::new(&fem, &config)?;
    let mut output = FemOutput::new(&fem, None, None, None)?;

    // solution
    let mut solver = FemSolverImplicit::new(&fem, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00, -5.310671749739340e-06,
        -4.211153193865630e-07, -5.041231308584792e-06,
        -7.221882229919595e-07, -3.342856970059000e-06,
        -4.211153193865630e-07, -1.644482664256252e-06,
         0.000000000000000e+00, -1.375042190378655e-06,
         0.000000000000000e+00, -4.288483503711751e-06,
         3.774202249726477e-07, -2.785714138358060e-06,
         0.000000000000000e+00, -1.282944773004366e-06,
         0.000000000000000e+00, -3.243122913926113e-06,
         2.708147610339125e-07, -2.873165860694200e-06,
         3.669995660268174e-07, -2.228571333618901e-06,
         2.708147610339123e-07, -1.583976767620262e-06,
         0.000000000000000e+00, -1.214019753311685e-06,
         0.000000000000000e+00, -2.217378283206796e-06,
         2.996313570447692e-07, -1.671428489956403e-06,
         0.000000000000000e+00, -1.125478696706005e-06,
         0.000000000000000e+00, -1.378773528167206e-06,
         1.370485290061350e-07, -1.298706937331023e-06,
         1.912334843273510e-07, -1.114285662470972e-06,
         1.370485290061355e-07, -9.298643811646859e-07,
         0.000000000000000e+00, -8.497977967747343e-07,
         0.000000000000000e+00, -6.453528854650004e-07,
         1.121994376475524e-07, -5.571428285393078e-07,
         0.000000000000000e+00, -4.689327716136146e-07,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(&state.uu, uu_correct, 1e-12);
    Ok(())
}
