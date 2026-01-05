use gemlab::graph::GraphUnd;
use gemlab::prelude::*;
use plotpy::Plot;
use pmsim::fem::PostProcMemo;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;
use russell_lab::vec_inner;
use russell_lab::Vector;
use tritet::Trigen;

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_ell_shape_simple";
const RENUMBER_POINTS: bool = false;
const DRAW_TRIANGLES: bool = false;
const DRAW_MESH: bool = false;

#[test]
fn test_seep_ell_shape() -> Result<(), StrError> {
    // mesh
    let generate = false;
    let triangle = true;
    let finer = false;
    let finest = false;
    let o2 = true;
    let mesh = generate_or_read_mesh(generate, triangle, finer, finest, o2);

    // features
    let features = Features::new(&mesh, false);
    let inlet = features.search_edges(At::Y(2.0), |_| true)?;
    let outlet = features.search_edges(At::X(2.0), |_| true)?;

    // parameters
    let (kx, ky) = (1.0, 1.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&inlet, Dof::Phi, 25.0);
    essential.edges(&outlet, Dof::Phi, 0.0);

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_out_files(OUT_DIR, NAME, 1.0)
        .update_model_settings(1)
        .set_save_flux(true);

    // FEM state
    let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &schema, &config)?;

    // solution
    let mut solver = SolverOld::new(&mesh, &schema, &config, &essential, &natural)?;
    solver.solve_sys(&mut state, &mut results)?;

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // post-processing tool
    let (post, mut memo) = PostProc::new(OUT_DIR, NAME)?;

    // read last state
    let last = post.n_state() - 1;
    let state = post.read_state(last)?;

    // load mesh and find vertical section along the gap underneath the wall
    let mesh = post.mesh();
    let with_internal_edges = true; // need internal edges
    let cells_by_points = true; // use all cells surrounding a point for better extrapolation
    let features = Features::new(&mesh, with_internal_edges);
    let (_, top_cells, top_section) = features.search_features_2d(cells_by_points, At::Y(2.0), |_| true)?;
    let (_, mid_cells, mid_section) = features.search_features_2d(cells_by_points, At::X(1.0), |x| x[1] <= 1.0)?;
    let (_, rig_cells, rig_section) = features.search_features_2d(cells_by_points, At::X(2.0), |_| true)?;
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
    println!("Flux through top section = {}", q0);

    // flux through mid section
    let un1 = Vector::from(&[1.0, 0.0]);
    let q1 = calc_flux_through_flat_section(&post, &mut memo, &state, &mid_cells, &mid_section, &un1)?;
    println!("Flux through mid section = {}", q1);

    // flux through right section
    let un2 = Vector::from(&[1.0, 0.0]); // same as mid section
    let q2 = calc_flux_through_flat_section(&post, &mut memo, &state, &rig_cells, &rig_section, &un2)?;
    println!("Flux through rig section = {}", q2);

    // Results from the finest Quad/O2 mesh:
    // Flux through top section = -9.771429623018324
    // Flux through mid section = 9.490741439420841
    // Flux through rig section = 9.771429622853258

    approx_eq(q0, 9.77, 0.15);

    // write Paraview files
    let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    println!("Paraview File: {}", path_pvd);

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

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool, triangle: bool, finer: bool, finest: bool, o2: bool) -> Mesh {
    let mut key = if finer {
        if finest {
            "finest".to_string()
        } else {
            "finer".to_string()
        }
    } else {
        "coarse".to_string()
    };
    if o2 {
        key = format!("{}_o2", key);
    }
    if generate {
        let mut mesh = if triangle {
            // set perimeter
            let mut gen = Trigen::new(7, Some(7), Some(2), None).unwrap();
            gen.set_point(0, 0, 0.0, 0.0).unwrap();
            gen.set_point(1, 0, 2.0, 0.0).unwrap();
            gen.set_point(2, 0, 2.0, 1.0).unwrap();
            gen.set_point(3, 0, 1.0, 1.0).unwrap();
            gen.set_point(4, 0, 1.0, 2.0).unwrap();
            gen.set_point(5, 0, 0.0, 2.0).unwrap();
            gen.set_point(6, 0, 1.0, 0.0).unwrap();
            gen.set_segment(0, 0, 0, 1).unwrap();
            gen.set_segment(1, 1, 1, 2).unwrap();
            gen.set_segment(2, 0, 2, 3).unwrap();
            gen.set_segment(3, 0, 3, 4).unwrap();
            gen.set_segment(4, 4, 4, 5).unwrap();
            gen.set_segment(5, 0, 5, 0).unwrap();
            gen.set_segment(6, 0, 6, 3).unwrap();
            gen.set_region(0, 1, 0.1, 0.1, None).unwrap();
            gen.set_region(1, 1, 1.1, 0.1, None).unwrap();

            // generate mesh
            let max_area = if finest {
                0.001
            } else if finer {
                0.01
            } else {
                0.2
            };
            gen.generate_mesh(false, o2, true, Some(max_area), None).unwrap();

            if DRAW_TRIANGLES {
                let mut plot = Plot::new();
                gen.draw_triangles(&mut plot, true, false, false, false, None, None, None);
                plot.set_figure_size_points(800.0, 800.0)
                    .save(&format!("{}/mesh_{}_tri.svg", OUT_DIR, NAME))
                    .unwrap();
            }

            Unstructured::from_trigen(&gen)
        } else {
            // set the geometry kind
            let kind = if o2 { GeoKind::Qua8 } else { GeoKind::Qua4 };

            // allocate blocks
            let mut block1 = Block::new(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).unwrap();
            let mut block2 = Block::new(&[[0.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]]).unwrap();
            let mut block3 = Block::new(&[[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]]).unwrap();

            // subdivide blocks
            let n = if finest {
                200
            } else if finer {
                20
            } else {
                2
            };
            block1.set_ndiv(&[n, n]).unwrap();
            block2.set_ndiv(&[n, n]).unwrap();
            block3.set_ndiv(&[n, n]).unwrap();
            let mesh1 = block1.subdivide(kind).unwrap();
            let mesh2 = block2.subdivide(kind).unwrap();
            let mesh3 = block3.subdivide(kind).unwrap();

            // join meshes
            join_meshes(&[&mesh1, &mesh2, &mesh3]).unwrap()
        };

        if RENUMBER_POINTS {
            GraphUnd::renumber_mesh(&mut mesh, true).unwrap();
        }
        mesh.check_all().unwrap();

        if DRAW_MESH {
            let coarse = !finer && !finest;
            let show_ids = coarse;
            let mut draw = Draw::new();
            draw.show_point_ids(show_ids)
                .show_cell_ids(show_ids)
                .show_cell_marker(false)
                .set_size(800.0, 800.0)
                .all(&mesh, &format!("{}/mesh_{}_{}.svg", OUT_DIR, NAME, key))
                .unwrap();
        }

        // write mesh
        mesh.write(&format!("{}/{}_{}.msh", OUT_DIR, NAME, key)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}_{}.msh", NAME, key)).unwrap()
    }
}
