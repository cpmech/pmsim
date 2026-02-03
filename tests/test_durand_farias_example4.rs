use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::analytical::ElastPlaneStrainFlexibleFoot;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;
use std::fmt::Write;

const NAME: &str = "test_durand_farias_example4";
const KIND: GeoKind = GeoKind::Qua4;
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = false;

// constants
const B: f64 = 2.5; // half-width of the flexible footing
const W: f64 = 15.0 * B; // width of the domain
const H: f64 = 15.0 * B; // height of the domain
const QN: f64 = 200.0; // magnitude of the distributed loading
const E: f64 = 1000.0; // Young's modulus
const NU: f64 = 0.0; // Poisson's coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_durand_farias_example4() -> Result<(), StrError> {
    run_test(true)?;
    run_test(false)?;
    Ok(())
}

fn run_test(lmm: bool) -> Result<(), StrError> {
    // mesh
    let mesh = generate_or_read_mesh(1, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(W), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let footing = features.search_edges(At::Y(H), |x| x[0] <= B)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic { young: E, poisson: NU },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.edges(&footing, Nbc::Qn, -QN);

    // configuration
    let key = if lmm { "_lmm" } else { "_sps" };
    let name = &format!("{}{}", NAME, key);
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(lmm)
        .set_out_files("/tmp/pmsim", name, 1.0);

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;

    //
    // data analysis -------------------------------------------------------------
    //

    // results
    let (post, mut memo) = PostProc::new("/tmp/pmsim", name)?;

    // features
    let left = features.search_edges(At::X(0.0), any_x)?;
    let left_cells = features.get_cells_via_2d_edges(&left);
    let (min, max) = mesh.get_cell_bounding_box(mesh.cells[left_cells[0]].id);
    let hdx = (max[0] - min[0]) / 2.0;

    // stresses
    let state = post.read_state(post.nstate() - 1)?;
    let gauss = post.gauss_stresses_patch(&mut memo, &state, &left_cells, |x, _, _| x < hdx)?;
    let nodal = post.nodal_stresses_patch(&mut memo, &state, &left_cells, |x, _, _| x < hdx)?;

    // verification
    let ana = ElastPlaneStrainFlexibleFoot {
        bb: B,
        hh: H,
        ww: W,
        qn: QN,
        young: E,
        poisson: NU,
    };
    let thin_line = format!("{:─^1$}", "", 6 * 8 + 4 * 3 + 5 * 2);
    let mut buf = String::new();
    writeln!(&mut buf, "\n{}", thin_line).unwrap();
    for i in 0..nodal.xx.len() {
        let x = nodal.xx[i];
        let y = nodal.yy[i];
        writeln!(
            &mut buf,
            "{:8.3} =? {:8.3}  │  {:8.3} =? {:8.3}  │  {:8.3} =? {:8.3}",
            nodal.txx[i],
            ana.stress(x, y).get(0, 0),
            nodal.tyy[i],
            ana.stress(x, y).get(1, 1),
            nodal.txy[i],
            ana.stress(x, y).get(0, 1)
        )
        .unwrap();
        approx_eq(f64::abs(nodal.txx[i] - ana.stress(x, y).get(0, 0)) / QN, 0.0, 0.25);
        approx_eq(f64::abs(nodal.tyy[i] - ana.stress(x, y).get(1, 1)) / QN, 0.0, 0.23);
        approx_eq(f64::abs(nodal.txy[i] - ana.stress(x, y).get(0, 1)) / QN, 0.0, 0.09);
    }
    writeln!(&mut buf, "{}", thin_line).unwrap();
    println!("{}", buf);
    let correct = r#"
──────────────────────────────────────────────────────────────────────
   5.355 =?   -0.025  │   -23.442 =?  -16.926  │    -0.000 =?   -0.000
   5.616 =?   -0.042  │   -24.391 =?  -20.152  │     0.085 =?   -0.000
   6.304 =?   -0.080  │   -27.534 =?  -24.887  │     0.177 =?   -0.000
   7.640 =?   -0.179  │   -33.648 =?  -32.498  │     0.408 =?   -0.000
   6.897 =?   -0.537  │   -44.297 =?  -46.662  │     0.577 =?   -0.000
  14.614 =?   -2.993  │   -76.212 =?  -81.117  │     0.915 =?   -0.000
 -46.706 =?  -90.037  │  -147.282 =? -191.896  │    15.255 =?   -0.000
-150.447 =? -200.000  │  -193.167 =? -200.000  │    17.513 =?   -0.000
──────────────────────────────────────────────────────────────────────
"#;
    assert_eq!(buf, correct);

    // figure
    if SAVE_FIGURE {
        let mut plot = Plot::new();
        let mut curve_ana = Curve::new();
        curve_ana
            .set_label("analytical (semi-infinite domain)")
            .set_line_style("--");
        let (ss, ll) = ana.get_normalized_syy_along_center(101);
        curve_ana.draw(&ss, &ll);

        let mut curve_gauss = Curve::new();
        curve_gauss
            .set_label("integration points")
            .set_marker_style("o")
            .set_marker_void(true)
            .set_line_style("None");
        let ll: Vec<_> = gauss.yy.iter().map(|y| y / B).collect();
        let ss: Vec<_> = gauss.tyy.iter().map(|sy| -sy / QN).collect();
        curve_gauss.draw(&ss, &ll);

        let mut curve_nodal = Curve::new();
        curve_nodal
            .set_label("nodal points (averaged)")
            .set_marker_style("s")
            .set_marker_size(5.0)
            .set_line_style("None");
        let ll: Vec<_> = nodal.yy.iter().map(|y| y / B).collect();
        let ss: Vec<_> = nodal.tyy.iter().map(|sy| -sy / QN).collect();
        curve_nodal.draw(&ss, &ll);

        plot.set_title("Durand and Farias Fig 15")
            .add(&curve_ana)
            .add(&curve_gauss)
            .add(&curve_nodal)
            .grid_labels_legend("Normalized stress: $-\\sigma_v/q_n$", "Normalized length: $y/B$");

        plot.save(&format!("/tmp/pmsim/{}_{}.svg", name, KIND.to_string()))?;
    }
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(att: i32, generate: bool) -> Mesh {
    let k_str = KIND.to_string();

    if generate {
        // generate mesh
        let (xa, xb, xc, na, nb) = (0.0, B, W, 4, 3);
        let (ya, yb, yc, ma, mb) = (0.0, H - B / 2.0, H, 6, 1);
        let mesh = Structured::rectangle(
            xa,
            Some(xb),
            xc,
            na,
            nb,
            &[ya, yb, yc],
            &[ma, mb],
            &[att, att],
            KIND,
            true,
        )
        .unwrap();
        mesh.check_all().unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(true)
            .show_cell_ids(true)
            .set_size(1000.0, 1000.0)
            .all(&mesh, &format!("/tmp/pmsim/mesh_{}_{}.svg", NAME, k_str))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}_{}.msh", NAME, k_str)).unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}_{}.msh", NAME, k_str)).unwrap()
    }
}
