#![allow(unused)]

use gemlab::graph::GraphUnd;
use gemlab::prelude::*;
use plotpy::Plot;
use pmsim::prelude::*;
use pmsim::StrError;
use tritet::Trigen;

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_craig_2d9";
const GENERATE_MESH: bool = true;

// geometry constants
const G: f64 = 2.25; // gap between the bedrock and pile wall
const P: f64 = 6.0; // pile wall height
const L: f64 = G + P + 2.0; // left height
const R: f64 = G + P; // right height
const M: f64 = 1.0 * L; // left distance (it should be > 1.5 L)
const N: f64 = 5.5 / 2.0; // right distance
const T: f64 = 0.5; // thickness of sheet pile wall

// number of divisions (coarse)
const NYM_COARSE: usize = 2; // number of divisions in M
const NYN_COARSE: usize = 1; // number of divisions in N
const NYL_COARSE: usize = 1; // number of divisions for the extra height on the left side
const NYG_COARSE: usize = 1; // number of divisions for the gap
const NYP_COARSE: usize = 2; // number of divisions for the pile

// number of divisions (finer)
const NYM_FINER: usize = 4; // number of divisions in M
const NYN_FINER: usize = 3; // number of divisions in N
const NYL_FINER: usize = 2; // number of divisions for the extra height on the left side
const NYG_FINER: usize = 4; // number of divisions for the gap
const NYP_FINER: usize = 4; // number of divisions for the pile

#[test]
fn test_seep_craig_2d9() -> Result<(), StrError> {
    // mesh
    let finer = false;
    let o2 = false;
    // let mesh = generate_or_read_mesh(GENERATE_MESH, finer, o2);
    let mesh = generate_or_read_mesh_simple();
    // let mesh = generate_or_read_mesh_simple_tri();

    // features
    let features = Features::new(&mesh, false);
    // let left_top = features.search_edges(At::Y(L), |x| x[0] <= M)?;
    // let right_top = features.search_edges(At::Y(R), |x| x[0] >= M + T)?;
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

    /*
    let gap_section = features.search_edges(At::X(M), |x| x[1] <= G)?;
    let cells_around_gap_section = features.get_cells_via_2d_edges(&gap_section);
    let gap_cells: Vec<_> = cells_around_gap_section
        .iter()
        .filter(|&&c| {
            for a in &mesh.cells[c].points {
                if mesh.points[*a].coords[0] < M {
                    return false;
                }
            }
            true
        })
        .copied()
        .collect();
    println!("section = {}", gap_section);
    println!("gap_cells = {:?}", gap_cells);
    */

    // analysis
    let state = post.read_state(post.n_state() - 1)?;
    let gauss = post.gauss_fluxes_patch(&mut memo, &state, &[6, 7], Dof::Phi, |_, _, _| true)?;
    // println!("wx = {:?}", gauss.vvx);
    // println!("wy = {:?}", gauss.vvy);

    let vv = post.nodal_fluxes_patch(
        &mut memo,
        &state,
        &[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
        Dof::Phi,
        |_, _, _| true,
    )?;
    // println!();
    // println!("wx = {:?}", vv.vvx);
    // println!("wy = {:?}", vv.vvy);
    // println!();
    for p in 0..mesh.points.len() {
        println!("point {:>2}: vx = {:15.10}, vy = {:15.10}", p, vv.vx(p), vv.vy(p));
    }

    // write Paraview files
    let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    println!("Paraview File: {}", path_pvd);

    // done
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool, finer: bool, o2: bool) -> Mesh {
    if generate {
        let kind = if o2 { GeoKind::Qua8 } else { GeoKind::Qua4 };

        // set blocks on the left-hand side
        let mut block1 = Block::new(&[[0.0, 0.0], [M, 0.0], [M, G], [0.0, G]]).unwrap();
        let mut block2 = Block::new(&[[0.0, G], [M, G], [M, G + P], [0.0, G + P]]).unwrap();
        let mut block3 = Block::new(&[[0.0, G + P], [M, G + P], [M, L], [0.0, L]]).unwrap();

        // set block under the sheet pile wall
        let mut block4 = Block::new(&[[M, 0.0], [M + T, 0.0], [M + T, G], [M, G]]).unwrap();

        // set blocks on the right-hand side
        let mut block5 = Block::new(&[[M + T, 0.0], [M + T + N, 0.0], [M + T + N, G], [M + T, G]]).unwrap();
        let mut block6 = Block::new(&[[M + T, G], [M + T + N, G], [M + T + N, R], [M + T, R]]).unwrap();

        // set number of division
        if finer {
            block1.set_ndiv(&[NYM_FINER, NYG_FINER]).unwrap();
            block2.set_ndiv(&[NYM_FINER, NYP_FINER]).unwrap();
            block3.set_ndiv(&[NYM_FINER, NYL_FINER]).unwrap();
            block4.set_ndiv(&[1, NYG_FINER]).unwrap();
            block5.set_ndiv(&[NYN_FINER, NYG_FINER]).unwrap();
            block6.set_ndiv(&[NYN_FINER, NYP_FINER]).unwrap();
        } else {
            block1.set_ndiv(&[NYM_COARSE, NYG_COARSE]).unwrap();
            block2.set_ndiv(&[NYM_COARSE, NYP_COARSE]).unwrap();
            block3.set_ndiv(&[NYM_COARSE, NYL_COARSE]).unwrap();
            block4.set_ndiv(&[1, NYG_COARSE]).unwrap();
            block5.set_ndiv(&[NYN_COARSE, NYG_COARSE]).unwrap();
            block6.set_ndiv(&[NYN_COARSE, NYP_COARSE]).unwrap();
        }

        // generate mesh
        let mesh1 = block1.subdivide(kind).unwrap();
        let mesh2 = block2.subdivide(kind).unwrap();
        let mesh3 = block3.subdivide(kind).unwrap();
        let mesh4 = block4.subdivide(kind).unwrap();
        let mesh5 = block5.subdivide(kind).unwrap();
        let mesh6 = block6.subdivide(kind).unwrap();
        let mut mesh = join_meshes(&[&mesh1, &mesh2, &mesh3, &mesh4, &mesh5, &mesh6]).unwrap();
        GraphUnd::renumber_mesh(&mut mesh, true).unwrap();
        mesh.check_all().unwrap();

        // draw figure
        let mut fig = Figure::new();
        fig.show_point_ids(true)
            .show_cell_ids(true)
            .show_cell_att(false)
            .size(800.0, 800.0)
            .draw(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("{}/{}.msh", OUT_DIR, NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}

fn generate_or_read_mesh_simple() -> Mesh {
    let mut block1 = Block::new(&[[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]).unwrap();
    let mut block2 = Block::new(&[[0.0, 1.0], [1.0, 1.0], [1.0, 2.0], [0.0, 2.0]]).unwrap();
    let mut block3 = Block::new(&[[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]]).unwrap();
    block1.set_ndiv(&[2, 2]).unwrap();
    block2.set_ndiv(&[2, 2]).unwrap();
    block3.set_ndiv(&[2, 2]).unwrap();
    let mesh1 = block1.subdivide(GeoKind::Qua4).unwrap();
    let mesh2 = block2.subdivide(GeoKind::Qua4).unwrap();
    let mesh3 = block3.subdivide(GeoKind::Qua4).unwrap();
    let mesh = join_meshes(&[&mesh1, &mesh2, &mesh3]).unwrap();
    mesh.check_all().unwrap();
    Figure::new()
        .show_point_ids(true)
        .show_cell_ids(true)
        .show_cell_att(false)
        .size(800.0, 800.0)
        .draw(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))
        .unwrap();
    mesh
}

fn generate_or_read_mesh_simple_tri() -> Mesh {
    let mut gen = Trigen::new(6, Some(6), Some(1), None).unwrap();
    gen.set_point(0, 0, 0.0, 0.0).unwrap();
    gen.set_point(1, 0, 2.0, 0.0).unwrap();
    gen.set_point(2, 0, 2.0, 1.0).unwrap();
    gen.set_point(3, 0, 1.0, 1.0).unwrap();
    gen.set_point(4, 0, 1.0, 2.0).unwrap();
    gen.set_point(5, 0, 0.0, 2.0).unwrap();
    gen.set_segment(0, 0, 0, 1).unwrap();
    gen.set_segment(1, 1, 1, 2).unwrap();
    gen.set_segment(2, 0, 2, 3).unwrap();
    gen.set_segment(3, 0, 3, 4).unwrap();
    gen.set_segment(4, 4, 4, 5).unwrap();
    gen.set_segment(5, 0, 5, 0).unwrap();
    gen.set_region(0, 1, 0.1, 0.1, None).unwrap();
    gen.generate_mesh(false, false, true, Some(0.01), None).unwrap();

    let mut plot = Plot::new();
    // gen.draw_triangles(&mut plot, true, false, false, false, None, None, None);
    plot.set_figure_size_points(800.0, 800.0)
        .save(&format!("{}/mesh_{}_tri.svg", OUT_DIR, NAME))
        .unwrap();

    let mesh = Unstructured::from_trigen(&gen);
    mesh.check_all().unwrap();

    /*
    Figure::new()
        .show_point_ids(true)
        .show_cell_ids(true)
        .show_cell_att(false)
        .size(800.0, 800.0)
        .draw(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))
        .unwrap();
    */

    mesh
}
