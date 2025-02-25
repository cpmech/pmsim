use plotpy::Curve;
use plotpy::Plot;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{approx_eq, read_data};

// Arc-length control (ALC) for a 12-member truss in 3D
//
// TEST GOAL
//
// This test verifies the arc-length (AL) implementation in SolverImplicit
//
// MESH
//
// The mesh consists of a 12-member truss in 2D as Ex 3.2 of Ref #1.
//
// BOUNDARY CONDITIONS
//
// Fully fixed @ the "base"
// Downward forces @ the "top"
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1.0 for all members
// Cross-sectional area = 1.0 for all members
// Geometrically-nonlinear rods with Green-Lagrange strain measure
//
// REFERENCES
//
// * Kadapa C (2021) A simple extrapolated predictor for overcoming the starting and tracking
//   issues in the arc-length method for nonlinear structural mechanics,
//   Engineering Structures, 234:111755

const NAME: &str = "test_alc_3d_truss_12member";
const SAVE_FIGURE: bool = false;

#[test]
fn test_alc_3d_truss_12member() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::truss_12member_3d();

    // parameters
    let p1 = ParamRod {
        gnl: Some(GnlStrain::Green),
        density: 1.0,
        young: 1.0,
        area: 1.0,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();

    // essential boundary conditions
    let mut essential = Essential::new();
    for p in [0, 1, 2, 6, 7, 8] {
        essential
            .point(p, Dof::Ux, 0.0)
            .point(p, Dof::Uy, 0.0)
            .point(p, Dof::Uz, 0.0);
    }

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.point(3, Pbc::Fz, -1.5);
    natural.point(4, Pbc::Fz, -1.0);
    natural.point(5, Pbc::Fz, -1.5);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_incremental(101) // 100 (as in ref #1) + 1 (initial state)
        .set_arc_length_method(true)
        .set_ini_trial_load_factor(0.025)
        .set_tol_rr_abs(1e-6)
        .set_tol_mdu_rel(1e-12);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

    // FileIo for writing results
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solver
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural).unwrap();
    solver.solve(&mut state, &mut file_io).unwrap();

    // check the total number of converged iterations
    assert_eq!(solver.n_converged_iterations(), 100);

    // analyze results
    analyze_results()
}

fn analyze_results() -> Result<(), StrError> {
    // load results
    let (post, _) = PostProc::new("/tmp/pmsim", NAME)?;
    let n = post.n_state();
    let mut arr_ell = Vec::with_capacity(n);
    let mut arr_u = Vec::with_capacity(n);
    let mut arr_v = Vec::with_capacity(n);
    let mut arr_w = Vec::with_capacity(n);
    let equ = post.eq(3, Dof::Ux)?;
    let eqv = post.eq(3, Dof::Uz)?;
    let eqw = post.eq(4, Dof::Uz)?;
    for i in 0..n {
        let state = post.read_state(i)?;
        arr_ell.push(state.ell);
        arr_u.push(state.u[equ]);
        arr_v.push(state.u[eqv]);
        arr_w.push(state.u[eqw]);
    }

    // load reference results from Reference #1
    let reference = read_data(
        "data/arc_length/kadapa-truss-3d-12members-model2.txt",
        &["ell", "u", "v", "w"],
    )?;

    // check number of load factors
    assert_eq!(arr_ell.len(), reference["ell"].len());

    // check results
    for i in 0..reference["ell"].len() {
        approx_eq(arr_ell[i], reference["ell"][i], 1e-5);
        approx_eq(arr_u[i], reference["u"][i], 1e-5);
        approx_eq(arr_v[i], reference["v"][i], 1e-5);
        approx_eq(arr_w[i], reference["w"][i], 1e-5);
    }

    // plot results
    if SAVE_FIGURE {
        let mut curve1 = Curve::new();
        let mut curve2 = Curve::new();
        let mut curve3 = Curve::new();
        let mut curve1_ref = Curve::new();
        let mut curve2_ref = Curve::new();
        let mut curve3_ref = Curve::new();
        curve1
            .set_label("pmsim: u")
            .set_line_style("None")
            .set_marker_style("o")
            .set_marker_void(true);
        curve1_ref.set_label("Kadapa (2021)").set_line_style("-");
        curve2
            .set_label("pmsim: v")
            .set_line_style("None")
            .set_marker_style("o")
            .set_marker_void(true);
        curve2_ref.set_label("Kadapa (2021)").set_line_style("-");
        curve3
            .set_label("pmsim: w")
            .set_line_style("None")
            .set_marker_style("o")
            .set_marker_void(true);
        curve3_ref.set_label("Kadapa (2021)").set_line_style("-");
        curve1.draw(&arr_u, &arr_ell);
        curve2.draw(&arr_v, &arr_ell);
        curve3.draw(&arr_w, &arr_u);
        curve1_ref.draw(&reference["u"], &reference["ell"]);
        curve2_ref.draw(&reference["v"], &reference["ell"]);
        curve3_ref.draw(&reference["w"], &reference["u"]);
        let mut plot = Plot::new();
        plot.set_subplot(3, 1, 1)
            .add(&curve1_ref)
            .add(&curve1)
            .grid_labels_legend("u: |ux| displacement @ 3", "load factor")
            .set_subplot(3, 1, 2)
            .add(&curve2_ref)
            .add(&curve2)
            .grid_labels_legend("v: |uz| displacement @ 3", "load factor")
            .set_subplot(3, 1, 3)
            .add(&curve3_ref)
            .add(&curve3)
            .grid_labels_legend("w: uz displacement @ 4", "v: uz displacement @ 3")
            .set_figure_size_points(600.0, 900.0)
            .save(&format!("/tmp/pmsim/{}.svg", NAME))
            .unwrap();
    }
    Ok(())
}
