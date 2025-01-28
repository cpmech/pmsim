use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::analytical::FlexibleFooting2d;
use pmsim::prelude::*;
use russell_lab::*;
use std::collections::HashMap;

const NAME: &str = "test_durand_farias_example4";
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
    };
    let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))])?;

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
    let mut config = Config::new(&mesh);
    config.set_ngauss(att, NGAUSS);

    // FEM state
    let mut state = FemState::new(&fem, &config)?;
    let mut output = FileIo::new(&fem, None, None)?;

    // solution
    let mut solver = FemSolverImplicit::new(&fem, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // results
    let mut post = PostProc::new(&fem, &config);
    let mut gauss_x = Vec::new();
    let mut gauss_y = Vec::new();
    let mut gauss_syy = Vec::new();
    let mut map_nodal_count = HashMap::new();
    let mut map_nodal_syy = HashMap::new();
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
        let syy = post.stress(cell_id, &state, 1, 1)?;
        for p in 0..ngauss {
            if f64::abs(gcs[p][0] - x_min) < 1e-3 {
                gauss_x.push(gcs[p][0]);
                gauss_y.push(gcs[p][1]);
                gauss_syy.push(syy[p]);
            }
        }
        let syy = post.stress_nodal(cell_id, &state, 1, 1)?;
        let nnode = syy.dim();
        for m in 0..nnode {
            let nid = mesh.cells[cell_id].points[m];
            if mesh.points[nid].coords[0] < 1e-3 {
                map_nodal_count.entry(nid).and_modify(|v| *v += 1).or_insert(1_usize);
                map_nodal_syy.entry(nid).and_modify(|v| *v += syy[m]).or_insert(syy[m]);
            }
        }
    }
    let mut nodal_x = Vec::new();
    let mut nodal_y = Vec::new();
    let mut nodal_syy = Vec::new();
    for nid in map_nodal_syy.keys() {
        let count = map_nodal_count.get(&nid).unwrap();
        let correct = if *nid == 0 || *nid == 14 { 1 } else { 2 };
        assert_eq!(*count, correct);
        let syy = map_nodal_syy.get(nid).unwrap();
        nodal_x.push(mesh.points[*nid].coords[0]);
        nodal_y.push(mesh.points[*nid].coords[1]);
        nodal_syy.push(*syy / (*count as f64));
    }

    // verification
    let ana = FlexibleFooting2d {
        bb: B,
        hh: H,
        ww: W,
        qn: QN,
        young: E,
        poisson: NU,
    };
    for i in 0..nodal_x.len() {
        let x = nodal_x[i];
        let y = nodal_y[i];
        let syy_num = nodal_syy[i];
        let syy_ana = ana.stress(x, y).get(1, 1);
        let rel_err = 100.0 * f64::abs(syy_num - syy_ana) / f64::abs(syy_ana);
        // println!("{} =? {} => {}", syy_num, syy_ana, rel_err);
        assert!(rel_err < 39.0); // yes, high percentage because the analytical solution is for an infinite medium
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
        let ll: Vec<_> = nodal_y.iter().map(|y| y / B).collect();
        let ss: Vec<_> = nodal_syy.iter().map(|sy| -sy / QN).collect();
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
