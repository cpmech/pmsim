use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::analytical::FlexibleFooting2d;
use pmsim::prelude::*;
use russell_lab::*;

const NAME: &str = "test_durand_farias_example4";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = true;

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
    // mesh
    let att = 1;
    let kind = GeoKind::Qua4;
    let mesh = generate_or_read_mesh(att, kind, GENERATE_MESH);

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
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&footing, Nbc::Qn, -QN);

    // configuration
    let config = Config::new(&mesh);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // results
    let mut post = PostProc::new(&mesh, &base);
    let mut gauss_x = Vec::new();
    let mut gauss_y = Vec::new();
    let mut gauss_syy = Vec::new();
    for edge in &left.all {
        let attached_cells = features.get_cells_via_2d_edge(edge);
        let cell_id = attached_cells[0]; // only one cell because the edge is on boundary
        let gcs = post.gauss_coords(cell_id)?;
        let ngauss = gcs.len();
        let mut x_min = gcs[0][0];
        for p in 0..ngauss {
            if gcs[p][0] < x_min {
                x_min = gcs[p][0];
            }
        }
        let sig = post.gauss_stress(cell_id, &state)?;
        for p in 0..ngauss {
            if f64::abs(gcs[p][0] - x_min) < 1e-3 {
                gauss_x.push(gcs[p][0]);
                gauss_y.push(gcs[p][1]);
                gauss_syy.push(sig.get(p, 0));
            }
        }
    }

    // extrapolated results
    let cell_ids = features.get_cells_via_2d_edges(&left);
    let nodal = post.nodal_stresses(&cell_ids, &state, |_, x, _, _| x < 1e-3)?;

    // verification
    let ana = FlexibleFooting2d {
        bb: B,
        hh: H,
        ww: W,
        qn: QN,
        young: E,
        poisson: NU,
    };
    for i in 0..nodal.xx.len() {
        let x = nodal.xx[i];
        let y = nodal.yy[i];
        println!(
            "{:8.3} =? {:8.3}   #   {:8.3} =? {:8.3}   #   {:8.3} =? {:8.3}",
            nodal.txx[i],
            ana.stress(x, y).get(0, 0),
            nodal.tyy[i],
            ana.stress(x, y).get(1, 1),
            nodal.txy[i],
            ana.stress(x, y).get(0, 1)
        );
    }

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
        let ll: Vec<_> = gauss_y.iter().map(|y| y / B).collect();
        let ss: Vec<_> = gauss_syy.iter().map(|sy| -sy / QN).collect();
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

        plot.save(&format!("{}/{}_{}.svg", DEFAULT_TEST_DIR, NAME, kind.to_string()))?;
    }
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(att: usize, kind: GeoKind, generate: bool) -> Mesh {
    let k_str = kind.to_string();

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
            kind,
            true,
        )
        .unwrap();
        mesh.check_all().unwrap();

        // draw figure
        if SAVE_FIGURE {
            let mut opt = Figure::new();
            opt.with_point_marker = false;
            // opt.point_dots = true;
            opt.point_ids = true;
            opt.cell_ids = true;
            opt.with_cell_att = false;
            opt.figure_size = Some((1000.0, 1000.0));
            mesh.draw(
                Some(opt),
                &format!("{}/mesh_{}_{}.svg", DEFAULT_TEST_DIR, NAME, k_str),
                |_, _| {},
            )
            .unwrap();
        }

        // write mesh
        mesh.write_json(&format!("{}/{}_{}.json", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read_json(&format!("data/meshes/{}_{}.json", NAME, k_str)).unwrap()
    }
}
