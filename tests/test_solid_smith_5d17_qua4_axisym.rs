use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use russell_lab::*;

// Smith's Example 5.17 (Figure 5.17) on page 187
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies an axisymmetric equilibrium problem.
//
// MESH
//
//               1.0 kN/m²
//          ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
//   0.0    0------3----------6-------------------9
//       Ux | (0)  |   (2)    |        (4)        | Ux
//       F  | [1]  |   [1]    |        [1]        | F
//  -4.0 I  1------4----------7------------------10 I
//       X  | (1)  |   (3)    |        (5)        | X
//       E  | [2]  |   [2]    |        [2]        | E
// -10.0 D  2------5----------8------------------11 D
//         0.0    4.0       10.0                30.0
//                      Ux and Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix left edge horizontally
// Fix right edge horizontally
// Fix bottom edge horizontally and vertically
// Concentrated load (Fy) on points 0, 3, 6, equal to
// -2.6667, -23.3333, -24.0, respectively
// Distributed load Qn = -1 on top edge with x ≤ 4
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Upper layer: Young = 100, Poisson = 0.3
// Lower layer: Young = 1000, Poisson = 0.45
// Plane-strain
// NOTE: using 9 integration points

const VERBOSE_LEVEL: usize = 0;

#[test]
fn test_solid_smith_5d17_qua4_axisym() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d17_qua4();

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(30.0), any_x)?;
    let bottom = features.search_edges(At::Y(-10.0), any_x)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 100.0,
            poisson: 0.3,
        },
        ngauss: Some(9),
    };
    let p2 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 1000.0,
            poisson: 0.45,
        },
        ngauss: Some(9),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1)), (2, Elem::Solid(p2))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .points(&[0], Pbc::Fy, -2.6667)
        .points(&[3], Pbc::Fy, -23.3333)
        .points(&[6], Pbc::Fy, -24.0);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_alt_bb_matrix_method(true).set_axisymmetric();

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim/", "test_solid_smith_5d17_qua4_axisym")?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00, -3.176060471584131e-02,
        0.000000000000000e+00, -3.231272276712264e-03,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.394996413258417e-03, -3.990499561989564e-02,
        1.164812317196636e-03, -2.497995580677533e-03,
        0.000000000000000e+00,  0.000000000000000e+00,
        1.703556484195065e-03, -6.045921606453852e-03,
        1.330158368942836e-03, -4.421423576906877e-04,
        0.000000000000000e+00,  0.000000000000000e+00,
        0.000000000000000e+00,  2.587551628166839e-03,
        0.000000000000000e+00,  3.090608328409013e-04,
        0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(&state.uu, uu_correct, 1e-9);

    #[rustfmt::skip]
    let kk_e0_ref = Matrix::from(&[
        [ 4.145299153899789e+02, -6.410256675972500e+00,  1.282051134704921e+01,  6.410256914773158e+00, -1.282051246145228e+01, -3.205128258214119e+01,  1.880341878109254e+02,  3.205128234334054e+01],
        [-6.410256675972496e+00,  7.051282113389323e+01, -1.282051247519651e+01,  1.923076855330780e+01, -6.410256436770716e+01, -5.769230733324714e+01, -3.205128234334053e+01, -3.205128235395389e+01],
        [ 1.282051134704917e+01, -1.282051247519651e+01,  2.948717975954594e+02, -1.025641036258208e+02,  8.974358861192951e+01,  5.128205173331020e+01, -1.282051246145228e+01,  6.410256436770715e+01],
        [ 6.410256914773158e+00,  1.923076855330780e+01, -1.025641036258208e+02,  1.602564132190665e+02, -5.128205173331020e+01, -1.217948744391272e+02,  3.205128258214119e+01, -5.769230733324714e+01],
        [-1.282051246145228e+01, -6.410256436770716e+01,  8.974358861192951e+01, -5.128205173331020e+01,  2.948717975954594e+02,  1.025641036258208e+02,  1.282051134704917e+01,  1.282051247519651e+01],
        [-3.205128258214119e+01, -5.769230733324714e+01,  5.128205173331020e+01, -1.217948744391272e+02,  1.025641036258208e+02,  1.602564132190665e+02, -6.410256914773157e+00,  1.923076855330780e+01],
        [ 1.880341878109254e+02, -3.205128234334054e+01, -1.282051246145228e+01,  3.205128258214120e+01,  1.282051134704920e+01, -6.410256914773157e+00,  4.145299153899790e+02,  6.410256675972501e+00],
        [ 3.205128234334053e+01, -3.205128235395389e+01,  6.410256436770715e+01, -5.769230733324714e+01,  1.282051247519651e+01,  1.923076855330780e+01,  6.410256675972497e+00,  7.051282113389323e+01],
    ]);
    let e0 = &solver.elements.all[0];
    mat_approx_eq(&e0.kke, &kk_e0_ref, 1e-5);

    // compare the results with the reference
    let tol_displacement = 7e-10;
    let tol_stress = 8e-8;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SGM,
        "data/sgm/sgm_5d17_ref.json",
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);
    Ok(())
}
