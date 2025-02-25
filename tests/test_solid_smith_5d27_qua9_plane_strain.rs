use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::vec_approx_eq;

// Smith's Example 5.27 (Figure 5.27) on page 200
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a plane-strain simulation with Qua9 elements and full integration.
// NOTE: This Example is similar to Example 5.15, with the difference being Qua9 elements.
//
// MESH
//
//           1.0 kN/m²
//          ↓↓↓↓↓↓↓↓↓↓↓
//  0.0     0----1----2----3----4
//          |         |         |
//          5    6    7    8    9
//          |         |         |
// -3.0 Ux 10---11---12---13---14 Ux
//      F   |         |         | F
//      I  15   16   17   18   19 I
//      X   |         |         | X
// -6.0 E  20---21---22---23---24 E
//      D   |         |         | D
//         25   26   27   28   29
//          |         |         |
// -9.0    30---31---32---33---34
//         0.0       3.0       6.0
//             Ux and Uy FIXED
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

#[test]
fn test_solid_smith_5d27_qua9_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d27_qua9();

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
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

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
    let config = Config::new(&mesh);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00, -5.298917281673461e-06,
        -4.003799238378697e-07, -4.988343783990114e-06,
        -6.190173499951097e-07, -3.342857278217464e-06,
        -4.003799238378709e-07, -1.697370595362790e-06,
         0.000000000000000e+00, -1.386797274761462e-06,
         0.000000000000000e+00, -4.307432451606995e-06,
         1.856221968116303e-07, -3.910784322389835e-06,
         3.166914493685104e-07, -2.785714335355316e-06,
         1.856221968116306e-07, -1.660644261280014e-06,
         0.000000000000000e+00, -1.263996219103632e-06,
         0.000000000000000e+00, -3.169857671727359e-06,
         2.748430716746047e-07, -2.899079293310467e-06,
         3.776886772815824e-07, -2.228571506265132e-06,
         2.748430716746045e-07, -1.558063638804644e-06,
         0.000000000000000e+00, -1.287285340802904e-06,
         0.000000000000000e+00, -2.195216731858149e-06,
         2.109315386008372e-07, -2.046774420503323e-06,
         3.051479707607178e-07, -1.671428599981276e-06,
         2.109315386008371e-07, -1.296082730930499e-06,
         0.000000000000000e+00, -1.147640468104400e-06,
         0.000000000000000e+00, -1.373100658874431e-06,
         1.368995875875396e-07, -1.299314030835482e-06,
         1.947637253994351e-07, -1.114285752721248e-06,
         1.368995875875399e-07, -9.292574356333917e-07,
         0.000000000000000e+00, -8.554708465680631e-07,
         0.000000000000000e+00, -6.492347084869729e-07,
         7.643251179936223e-08, -6.245070433474793e-07,
         1.104696707838937e-07, -5.571428666315606e-07,
         7.643251179936240e-08, -4.897786738259930e-07,
         0.000000000000000e+00, -4.650510247761482e-07,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(&state.u, uu_correct, 3e-13);
    Ok(())
}
