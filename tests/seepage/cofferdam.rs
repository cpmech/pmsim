use gemlab::mesh::Blocks2d;
use gemlab::prelude::*;
use pmsim::fem::PostProcMemo;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{approx_eq, vec_inner, Vector};

const OUT_DIR: &str = "/tmp/pmsim/seepage";
const NAME: &str = "cofferdam";
const SAVE_FIGURE: bool = false;

fn draw_mesh(mesh: &Mesh, input: bool, tri: bool) -> Result<(), StrError> {
    let x_ticks = [-13.6, -4.6, 0.0, 4.4, 13.6];
    let x_labels = x_ticks.iter().map(|x| format!("{}", x)).collect::<Vec<_>>();
    let y_ticks = [0.0, 5.0, 7.5, 10.0];
    let y_labels = y_ticks.iter().map(|y| format!("{}", y)).collect::<Vec<_>>();
    let mut draw = Draw::new();
    draw.show_edge_markers(true)
        // .show_point_ids(true)
        // .show_cell_ids(true)
        .set_m_normal_vector_marker(0.025)
        .set_size(1200.0, 600.0)
        .set_view_flag(false)
        .extra(|plot, before| {
            if !before {
                plot.set_ticks_x_labels(&x_ticks, &x_labels)
                    .set_ticks_y_labels(&y_ticks, &y_labels)
                    .set_range(-15.0, 15.0, -1.0, 11.0);
            }
        });
    let path = if input {
        &format!("{}/mesh_in_{}.svg", OUT_DIR, NAME)
    } else {
        let key = if tri { "tri" } else { "qua" };
        &format!("{}/mesh_{}_{}.svg", OUT_DIR, NAME, key)
    };
    draw.all(&mesh, path)
}

#[test]
fn cofferdam() -> Result<(), StrError> {
    let use_input_mesh = false;
    let triangles = false;
    let o2_triangles = true;
    // mesh
    let mesh_in = Mesh::read("data/meshes/cofferdam.msh")?;
    if SAVE_FIGURE {
        draw_mesh(&mesh_in, true, false)?;
    }
    let mut blocks2d = Blocks2d {
        points: mesh_in.points.iter().map(|p| (p.coords[0], p.coords[1])).collect(),
        regions: mesh_in
            .cells
            .iter()
            .map(|c| (c.marker, c.points[0], c.points[1], c.points[2], c.points[3]))
            .collect(),
        div_weights: mesh_in.cells.iter().map(|_| (vec![1.0], vec![1.0])).collect(),
        edge_constraints: mesh_in.cells.iter().map(|_| None).collect(),
        marked_edges: mesh_in.marked_edges.clone(),
    };
    let mesh = if use_input_mesh {
        mesh_in
    } else if triangles {
        let holes = Vec::new();
        let global_max_area = Some(0.1);
        // with global_max_area = 0.001:
        // ndof       = 772075 │ dim(K)     = (773236,773236) │ genie = Umfpack
        // n_lagrange = 1161   │ nnz_sup(K) = 13864914        │
        // neq_total  = 773236 │ sym(K)     = YesFull         │
        // Flux through top section = 9.899e-7
        // Flux through mid section = 9.796e-7
        // Flux through rig section = 9.787e-7
        // q1 + q2 = 1.958e-6
        let mesh = Unstructured::call_trigen(&mesh_in, &holes, o2_triangles, None, global_max_area, None, false)?;
        if SAVE_FIGURE {
            draw_mesh(&mesh, false, true)?;
        }
        mesh
    } else {
        blocks2d.div_weights[1].1 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[1].0 = vec![3.0, 2.0, 1.0, 1.0, 1.0];
        blocks2d.div_weights[4].0 = vec![3.0, 2.0, 1.0, 1.0, 1.0];
        blocks2d.div_weights[4].1 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[2].1 = vec![3.0, 2.0, 1.0, 1.0, 1.0];
        blocks2d.div_weights[2].0 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[3].1 = vec![3.0, 2.0, 1.0, 1.0, 1.0];
        blocks2d.div_weights[3].0 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[5].0 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[6].1 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[7].0 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[7].1 = vec![1.0, 1.0, 2.0];
        blocks2d.div_weights[8].1 = vec![1.0, 1.0, 1.0];
        blocks2d.div_weights[8].0 = vec![1.0, 1.0, 2.0];
        blocks2d.div_weights[0].1 = vec![1.0, 1.0, 2.0];
        blocks2d.div_weights[0].0 = vec![1.0, 1.0];
        blocks2d.div_weights[9].0 = vec![1.0, 1.0, 2.0];
        blocks2d.div_weights[9].1 = vec![1.0, 1.0];
        let mesh = Structured::from_blocks_2d(&blocks2d, GeoKind::Qua4, false)?;
        if SAVE_FIGURE {
            draw_mesh(&mesh, false, false)?;
        }
        mesh
    };

    // features
    let with_internal_edges = true; // for data analysis
    let features = Features::new(&mesh, with_internal_edges);
    let inlet = features.search_marked_edges(20);
    let outlet = features.search_marked_edges(10);

    // parameters
    let (kx, ky) = (4e-7, 4e-7);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&inlet, Dof::Phi, 13.0);
    ebc.edges(&outlet, Dof::Phi, 7.5);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .lagrange_mult_method(true)
        .out_files(OUT_DIR, NAME)
        .update_model_settings(1)
        .set_save_flux(true);

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;

    //
    // data analysis -------------------------------------------------------------
    //

    // post-processing tool
    let (post, mut memo) = PostProc::new(OUT_DIR, NAME)?;

    // read last state
    let last = post.nfile() - 1;
    let state = post.read_file(last)?;

    // load mesh and find vertical section along the gap underneath the wall
    let cells_by_points = true; // use all cells surrounding a point for better extrapolation
    let (_, top_cells, top_section) = features.search_features_2d(cells_by_points, At::Y(10.0), |x| x[0] <= -4.6)?;
    let (_, mid_cells, mid_section) = features.search_features_2d(cells_by_points, At::X(-4.6), |x| x[1] <= 5.0)?;
    let (_, rig_cells, rig_section) = features.search_features_2d(cells_by_points, At::X(4.6), |x| x[1] <= 5.0)?;
    if mesh.cells.len() < 30 {
        println!("top_cells = {:?}", top_cells);
        println!("mid_cells = {:?}", mid_cells);
        println!("rig_cells = {:?}", rig_cells);
        println!(
            "top_section = {:?}",
            top_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
        println!(
            "mid_section = {:?}",
            mid_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
        println!(
            "rig_section = {:?}",
            rig_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
    }

    // flux through top section
    let un0 = Vector::from(&[0.0, -1.0]);
    let q0 = calc_flux_through_flat_section(&post, &mut memo, &state, &top_cells, &top_section, &un0)?;
    println!("Flux through top section: q0 = {:.3e}", q0);

    // flux through mid section
    let un1 = Vector::from(&[1.0, 0.0]);
    let q1 = calc_flux_through_flat_section(&post, &mut memo, &state, &mid_cells, &mid_section, &un1)?;
    println!("Flux through mid section: q1 = {:.3e}", q1);

    // flux through right section
    let un2 = Vector::from(&[-1.0, 0.0]); // to the left now
    let q2 = calc_flux_through_flat_section(&post, &mut memo, &state, &rig_cells, &rig_section, &un2)?;
    println!("Flux through rig section: q2 = {:.3e}", q2);

    // print the total flux (should equal 2 * q0)
    println!("q1 + q2 = {:.3e} ({:.3e})", q1 + q2, 2.0 * q0);

    // check
    let m = 1e6;
    approx_eq(q1 * m, q2 * m, 1e-13);
    approx_eq((q1 + q2) * m, 2.0 * q0 * m, 0.5);
    approx_eq(q1 * 1e7, 9.0, 1.0); // should be around 9.0
    approx_eq(q2 * 1e7, 9.0, 1.0); // should be around 9.0

    // write Paraview files
    // let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    // println!("Paraview File: {}", path_pvd);

    // done
    Ok(())
}

fn calc_flux_through_flat_section(
    post: &PostProc,
    memo: &mut PostProcMemo,
    state: &FemState,
    cells: &Vec<CellId>,
    section: &Edges,
    unit_normal: &Vector,
) -> Result<f64, StrError> {
    let patch = post.nodal_fluxes_patch(memo, &state, &cells, Dof::Phi, |_, _, _| true)?;
    let (points, coords, ww) = post.values_along_edges_vec(&patch, &section)?;
    let mut area = 0.0;
    let vertical_section = if unit_normal[0] == 0.0 {
        assert!(f64::abs(unit_normal[1]) == 1.0, "Unit normal must be unitary");
        false
    } else if unit_normal[1] == 0.0 {
        assert!(f64::abs(unit_normal[0]) == 1.0, "Unit normal must be unitary");
        true
    } else {
        panic!("Unit normal must be aligned with either x or y axis");
    };
    for i in 1..points.len() {
        let dx = coords[i][0] - coords[i - 1][0];
        let dy = coords[i][1] - coords[i - 1][1];
        if vertical_section {
            assert!(f64::abs(dx) < 1e-14, "Section must be vertical");
        } else {
            assert!(f64::abs(dy) < 1e-14, "Section must be horizontal");
        }
        let delta = if vertical_section { dy } else { dx };
        let w_dot_un_prev = vec_inner(&ww[i - 1], &unit_normal);
        let w_dot_un = vec_inner(&ww[i], &unit_normal);
        area += delta * (w_dot_un_prev + w_dot_un) / 2.0;
    }
    Ok(area)
}
