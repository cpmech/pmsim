use gemlab::graph::GraphUnd;
use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;
use russell_lab::vec_inner;
use russell_lab::Vector;

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_craig_2d9";

// Geometry:
//                      |←      5.5 m       →|
//
// ---WT--------------||                     ||
//           ↑        ||                     ||
//      A = 2.5 m     ||                     ||
//           ↓        ||                     ||
// ===================||                     ||===================
//           ↑        ||                     ||
//      D = 2.0 m     ||                     ||
//           ↓        ||                     ||
//          ---       ||===WT================||
//           ↑        ||                     ||
//           |        ||                     ||
//           |        ||                     ||
//      S = 6.0 m     ||                     ||
//           |        ||                     ||
//           |        ||                     ||
//           ↓        ||                     ||
//          ---       ||                     ||
//           ↑
//     F = 2.25 m
//           ↓
// ***********************IMPERMEABLE***DATUM*********************
//
// WT means water table
//
const A: f64 = 2.5; // height of water on the left side
const D: f64 = 2.0; // depth of the excavation
const S: f64 = 6.0; // pile wall height
const F: f64 = 2.25; // "flow" section; gap between the bedrock and pile wall
const L: f64 = D + S + F; // left height
const R: f64 = S + F; // right height
const M: f64 = 3.0 * L; // left distance
const N: f64 = 5.5 / 2.0; // right distance
const G: f64 = 0.01; // thickness of sheet pile wall

// number of divisions
const NYM: [usize; 3] = [2, 4, 8]; // number of divisions in M
const NYN: [usize; 3] = [1, 3, 6]; // number of divisions in N
const NYL: [usize; 3] = [1, 2, 4]; // number of divisions for the extra height on the left side
const NYF: [usize; 3] = [1, 4, 8]; // number of divisions for the "flow" section underneath the wall
const NYS: [usize; 3] = [2, 4, 8]; // number of divisions for the sheet pile wall

#[test]
fn test_seep_craig_2d9() -> Result<(), StrError> {
    // mesh
    let generate = true;
    let o2 = true;
    let index_refinement = 2; // 0, 1, or 2
    let mesh = generate_or_read_mesh(generate, o2, index_refinement);

    // features
    let features = Features::new(&mesh, false);
    let inlet = features.search_edges(At::Y(L), |x| x[0] <= M)?;
    let outlet = features.search_edges(At::Y(R), |x| x[0] >= M + G)?;

    // parameters
    let (kx, ky) = (2.6e-5, 2.6e-5); // m/s
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&inlet, Dof::Phi, L + A);
    essential.edges(&outlet, Dof::Phi, R);

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
    let mid_section = features.search_edges(At::X(M), |x| x[1] <= F)?;

    // read last state
    let last = post.n_state() - 1;
    let state = post.read_state(last)?;

    // extract fluxes along the mid section
    let mid_cell_ids = features.get_cells_via_2d_edges(&mid_section);
    let mid_patch = post.nodal_fluxes_patch(&mut memo, &state, &mid_cell_ids, Dof::Phi, |_, _, _| true)?;

    // approximate the integral of fluxes along the vertical section
    let (point_ids, coords, w_along_section) = post.values_along_edges_vec(&mid_patch, &mid_section)?;
    let mut area = 0.0;
    let un = Vector::from(&[1.0, 0.0]); // it's a vertical section
    for i in 1..point_ids.len() {
        let dy = coords[i][1] - coords[i - 1][1];
        let dot_prev = vec_inner(&w_along_section[i - 1], &un);
        let dot = vec_inner(&w_along_section[i], &un);
        area += dy * (dot_prev + dot) / 2.0;
    }
    let q_computed = 2.0 * area * 3600.0; // m³/h per unit length (multiply by 2 because of symmetry)
    let q_expected = 0.25; // m³/h as in the book
    println!("Flux through vertical middle section = {} ({})", q_computed, q_expected);
    approx_eq(q_computed, q_expected, 0.07);

    // write Paraview files
    let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    println!("\nParaview File: {}\n", path_pvd);
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool, o2: bool, index_refinement: usize) -> Mesh {
    if generate {
        assert!(index_refinement < 3, "index_refinement must be 0, 1, or 2");
        let kind = if o2 { GeoKind::Qua8 } else { GeoKind::Qua4 };
        let nym = NYM[index_refinement];
        let nyn = NYN[index_refinement];
        let nyl = NYL[index_refinement];
        let nyf = NYF[index_refinement];
        let nys = NYS[index_refinement];

        // set blocks on the left-hand side
        let mut block1 = Block::new(&[[0.0, 0.0], [M, 0.0], [M, F], [0.0, F]]).unwrap();
        let mut block2 = Block::new(&[[0.0, F], [M, F], [M, F + S], [0.0, F + S]]).unwrap();
        let mut block3 = Block::new(&[[0.0, F + S], [M, F + S], [M, L], [0.0, L]]).unwrap();

        // set block under the sheet pile wall
        let mut block4 = Block::new(&[[M, 0.0], [M + G, 0.0], [M + G, F], [M, F]]).unwrap();

        // set blocks on the right-hand side
        let mut block5 = Block::new(&[[M + G, 0.0], [M + G + N, 0.0], [M + G + N, F], [M + G, F]]).unwrap();
        let mut block6 = Block::new(&[[M + G, F], [M + G + N, F], [M + G + N, R], [M + G, R]]).unwrap();

        // set number of division
        block1.set_ndiv(&[nym, nyf]).unwrap();
        block2.set_ndiv(&[nym, nys]).unwrap();
        block3.set_ndiv(&[nym, nyl]).unwrap();
        block4.set_ndiv(&[1, nyf]).unwrap();
        block5.set_ndiv(&[nyn, nyf]).unwrap();
        block6.set_ndiv(&[nyn, nys]).unwrap();

        // subdivide blocks
        let mesh1 = block1.subdivide(kind).unwrap();
        let mesh2 = block2.subdivide(kind).unwrap();
        let mesh3 = block3.subdivide(kind).unwrap();
        let mesh4 = block4.subdivide(kind).unwrap();
        let mesh5 = block5.subdivide(kind).unwrap();
        let mesh6 = block6.subdivide(kind).unwrap();

        // join meshes
        let mut mesh = join_meshes(&[&mesh1, &mesh2, &mesh3, &mesh4, &mesh5, &mesh6]).unwrap();

        GraphUnd::renumber_mesh(&mut mesh, true).unwrap();
        mesh.check_all().unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(false)
            .show_cell_ids(false)
            .show_cell_marker(false)
            .set_view_flag(false)
            .set_size(1200.0, 800.0)
            .all(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("{}/{}.msh", OUT_DIR, NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
