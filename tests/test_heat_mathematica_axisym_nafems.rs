use gemlab::prelude::*;
use plotpy::Canvas;
use pmsim::{prelude::*, StrError};

// From Mathematica Heat Transfer Model Verification Tests
// (HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0002)
//
// NAFEMS benchmark test
//
// https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransferVerificationTests.html
//
// TEST GOAL
//
// This test verifies the steady heat equation with a localized flux boundary condition
//
// MESH (not-to-scale, not-equal-axis)
//
// 0.14     ||++++++++++++++++++++++++
//          ||   |    |    |    |    +
//          ||-----------------------+
//          ||   |    |    |    |    +
// 0.10  →→ |------------------------+  yb
//       →→ |    |    |    |    |    +
//       →→ |------------------------+
//       →→ |    |    |    |    |    +
//       →→ |------------------------+
//       →→ |    |    |    |    |    +
// 0.04  →→ |--------(#)-------------+  ya
//          ||   |    |    |    |    +
//          ||-----------------------+
//          ||   |    |    |    |    +
// 0.00     ||++++++++++++++++++++++++
//         0.02     0.04            0.10
//         rin      rref            rout
//
// '+' indicates sides with T = 273.15
// || means insulated
// →→ means flux with Qt = 5e5
// (#) indicates a reference point to check the results
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Temperature T = 273.15 on the top, bottom, and right edges
// Flux Qt = 5e5 on the middle-left edges from y=0.04 to y=0.10
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// No source
// Constant conductivity kx = ky = 52

const NAME: &str = "test_heat_mathematica_axisym_nafems";
const GENERATE_MESH: bool = false;
const REF_POINT_MARKER: PointMarker = -1;

#[test]
fn test_heat_mathematica_axisym_nafems() -> Result<(), StrError> {
    // geometry
    let (rin, rref, rout) = (0.02, 0.04, 0.1);
    let (ya, yb, h) = (0.04, 0.1, 0.14);

    // mesh
    let mesh = generate_or_read_mesh(rin, rref, rout, ya, yb, h, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let edges_temp = features.search_many_edges(&[At::Y(0.0), At::Y(h), At::X(rout)], any_x)?;
    let edges_flux = features.search_edges(At::X(rin), |x| x[1] >= ya && x[1] <= yb)?;

    // reference point
    let ref_point = mesh.search_first_marked_point(REF_POINT_MARKER, any_x)?;

    // reference solution
    let ref_temperature = 332.97;

    // parameters
    let (kx, ky) = (52.0, 52.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&edges_temp, Dof::T, 273.15);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&edges_flux, Nbc::Qt, 5e5);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_axisymmetric();

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check
    let eq = base.equations.eq(ref_point, Dof::T).unwrap();
    let rel_err = f64::abs(state.uu[eq] - ref_temperature) / ref_temperature;
    println!(
        "\nT = {:?}, reference = {:?}, rel_error = {:>.8} %",
        state.uu[eq],
        ref_temperature,
        rel_err * 100.0
    );
    assert!(rel_err < 1e-5);
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(rin: f64, rref: f64, rout: f64, ya: f64, yb: f64, h: f64, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let y = &[0.0, ya, yb, h];
        let att = &[1, 1, 1];
        let (na, nb, ny) = (4, 12, &[8, 12, 8]);
        let mut mesh = Structured::rectangle(rin, Some(rref), rout, na, nb, y, ny, att, GeoKind::Qua9, true).unwrap();

        // mark reference point
        let extract_all = true; // << needed to find interior point
        let features = Features::new(&mesh, extract_all);
        let ref_points = features.search_point_ids(At::XY(rref, ya), any_x).unwrap();
        assert_eq!(ref_points.len(), 1);
        mesh.points[ref_points[0]].marker = REF_POINT_MARKER;

        // write mesh
        mesh.write(&format!("{}/{}.msh", DEFAULT_TEST_DIR, NAME)).unwrap();

        // reference point
        let mut circle = Canvas::new();
        circle.set_face_color("red").set_edge_color("red");
        circle.draw_circle(rref, ya, 0.02 * (rout - rin));

        // configure plot
        let mut fig = Figure::new();
        fig.point_dots = true;
        fig.figure_size = Some((400.0, 600.0));
        fig.canvas_points.set_marker_size(3.0).set_marker_line_color("none");

        // generate figure
        mesh.draw(
            Some(fig),
            &format!("{}/{}.svg", DEFAULT_TEST_DIR, NAME),
            |plot, before| {
                if !before {
                    plot.add(&circle);
                }
            },
        )
        .unwrap();

        // return generated mesh
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
