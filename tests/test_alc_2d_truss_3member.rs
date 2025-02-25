use plotpy::Curve;
use plotpy::Plot;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{approx_eq, read_data};

// Arc-length control (ALC) for a 3-member truss in 2D
//
// TEST GOAL
//
// This test verifies the arc-length (AL) implementation in SolverImplicit
//
// MESH
//
// The mesh consists of a simple 3-member truss in 2D as Ex 3.1 of Ref #1.
//
//
// ```text
//       3
//       |        (#) indicates cell id
//   (2) | [2]    [#] indicates attribute id
//       |        The lines are ROD (Lin2) elements
//       1
// (0)  / \  (1)
// [1] /   \ [1]
//    /     \
//   0       2
// ```
//
// BOUNDARY CONDITIONS
//
// Fully fixed @ points 0 and 2
// Horizontally fixed @ point 3
// Downward force of -1.0 @ point 3
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1.0 for members with attribute id 1
// Young = 0.5 for members with attribute id 2
// Cross-sectional area = 1.0 for all members
// Geometrically-nonlinear rods with engineering strain measure
//
// REFERENCES
//
// * Kadapa C (2021) A simple extrapolated predictor for overcoming the starting and tracking
//   issues in the arc-length method for nonlinear structural mechanics,
//   Engineering Structures, 234:111755

const NAME: &str = "test_alc_2d_truss_3member";
const SAVE_FIGURE: bool = false;

#[test]
fn test_alc_2d_truss_3member() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::truss_3member_2d();

    // parameters
    let p1 = ParamRod {
        gnl: Some(GnlStrain::Eng),
        density: 1.0,
        young: 1.0,
        area: 1.0,
        ngauss: None,
    };
    let p2 = ParamRod {
        gnl: Some(GnlStrain::Eng),
        density: 1.0,
        young: 0.5,
        area: 1.0,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Rod(p1)), (2, Elem::Rod(p2))]).unwrap();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.point(0, Dof::Ux, 0.0);
    essential.point(0, Dof::Uy, 0.0);
    essential.point(2, Dof::Ux, 0.0);
    essential.point(2, Dof::Uy, 0.0);
    essential.point(3, Dof::Ux, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.point(3, Pbc::Fy, -1.0);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_incremental(51) // 50 (as in ref #1) + 1 (initial state)
        .set_arc_length_method(true)
        .set_ini_trial_load_factor(0.05)
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
    assert_eq!(solver.n_converged_iterations(), 46);

    // analyze results
    analyze_results()
}

fn analyze_results() -> Result<(), StrError> {
    // load results
    let (post, _) = PostProc::new("/tmp/pmsim", NAME)?;
    let n = post.n_state();
    let mut arr_ell = Vec::with_capacity(n);
    let mut arr_uy1 = Vec::with_capacity(n);
    let mut arr_uy3 = Vec::with_capacity(n);
    let eq1 = post.eq(1, Dof::Uy)?;
    let eq3 = post.eq(3, Dof::Uy)?;
    for i in 0..n {
        let state = post.read_state(i)?;
        arr_ell.push(state.ell);
        arr_uy1.push(state.u[eq1]);
        arr_uy3.push(state.u[eq3]);
    }

    // load reference results from Reference #1
    let reference = read_data(
        "data/arc_length/kadapa-truss-2d-3members-model1.txt",
        &["ell", "uy1", "uy3"],
    )?;

    // check number of load factors
    assert_eq!(arr_ell.len(), reference["ell"].len());

    // check results
    for i in 0..reference["ell"].len() {
        approx_eq(arr_ell[i], reference["ell"][i], 1e-5);
        approx_eq(arr_uy1[i], reference["uy1"][i], 1e-5);
        approx_eq(arr_uy3[i], reference["uy3"][i], 1e-5);
    }

    // plot results
    if SAVE_FIGURE {
        let mut curve1 = Curve::new();
        let mut curve2 = Curve::new();
        let mut curve1_ref = Curve::new();
        let mut curve2_ref = Curve::new();
        curve1
            .set_label("pmsim: uy1")
            .set_line_style("None")
            .set_marker_style("o")
            .set_marker_color("blue")
            .set_marker_line_color("blue");
        curve1_ref
            .set_label("Kadapa (2021)")
            .set_line_style("-")
            .set_line_color("blue");
        curve2
            .set_label("pmsim: uy3")
            .set_line_style("None")
            .set_marker_style("s")
            .set_marker_color("black")
            .set_marker_line_color("black");
        curve2_ref
            .set_label("Kadapa (2021)")
            .set_line_style("--")
            .set_line_color("black");
        curve1.draw(&arr_uy1, &arr_ell);
        curve2.draw(&arr_uy3, &arr_ell);
        curve1_ref.draw(&reference["uy1"], &reference["ell"]);
        curve2_ref.draw(&reference["uy3"], &reference["ell"]);
        let mut plot = Plot::new();
        plot.add(&curve1_ref)
            .add(&curve2_ref)
            .add(&curve1)
            .add(&curve2)
            .set_title("E(top element) = 0.5")
            .grid_labels_legend("vertical displacement", "load factor")
            .set_inv_x()
            .set_figure_size_points(600.0, 300.0)
            .save(&format!("/tmp/pmsim/{}.svg", NAME))
            .unwrap();
    }
    Ok(())
}
