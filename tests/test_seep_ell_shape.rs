use gemlab::graph::GraphUnd;
use gemlab::prelude::*;
use plotpy::Plot;
use pmsim::prelude::*;
use pmsim::StrError;
use tritet::Trigen;

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_ell_shape";
const RENUMBER_POINTS: bool = false;
const DRAW_TRIANGLES: bool = false;
const DRAW_MESH: bool = true;

#[test]
fn test_seep_ell_shape() -> Result<(), StrError> {
    // mesh
    let generate = true;
    let triangle = false;
    let finer = false;
    let finest = false;
    let o2 = false;
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
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&inlet, Dof::Phi, 25.0);
    essential.edges(&outlet, Dof::Phi, 0.0);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_out_files(OUT_DIR, NAME, 1.0)
        .update_model_settings(1)
        .set_save_flux(true);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &base, &config)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // post-processing tool
    let (post, mut memo) = PostProc::new(OUT_DIR, NAME)?;

    // load mesh and find vertical section along the gap underneath the wall
    let mesh = post.mesh();
    let with_internal_edges = true;
    let features = Features::new(&mesh, with_internal_edges); // need internal edges
    let mid_section = features.search_edges(At::X(1.0), |x| x[1] <= 1.0)?;
    println!("mid_section = {}", mid_section);
    let cell_ids = features.get_cells_via_2d_edges(&mid_section);
    let state = post.read_state(post.n_state() - 1)?;
    let q = post.integrate_flux_through_edges(&mut memo, &state, &mid_section, &cell_ids, Dof::Phi)?;
    println!("Flux through mid_section = {}", q);

    // write Paraview files
    let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    println!("Paraview File: {}", path_pvd);

    // done
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool, triangle: bool, finer: bool, finest: bool, o2: bool) -> Mesh {
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
            block1.set_ndiv(&[2, 2]).unwrap();
            block2.set_ndiv(&[2, 2]).unwrap();
            block3.set_ndiv(&[2, 2]).unwrap();
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
                .all(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))
                .unwrap();
        }

        // write mesh
        mesh.write(&format!("{}/{}.msh", OUT_DIR, NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
