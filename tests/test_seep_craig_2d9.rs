#![allow(unused)]

use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::Vector;

const NAME: &str = "test_seep_craig_2d9";
const GENERATE_MESH: bool = true;
const SAVE_FIGURE: bool = true;

#[test]
fn test_seep_craig_2d9() -> Result<(), StrError> {
    // mesh
    let mesh = generate_or_read_mesh(GENERATE_MESH);

    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool) -> Mesh {
    if generate {
        // constants
        const G: f64 = 2.25; // gap between the bedrock and pile wall
        const P: f64 = 6.0; // pile wall height
        const L: f64 = G + P + 2.0; // left height
        const R: f64 = G + P; // right height
        const M: f64 = 1.0 * L; // left distance (it should be > 1.5 L)
        const T: f64 = 0.001; // thickness of sheet pile wall
        const N: f64 = 5.5 / 2.0; // right distance
        let kind = GeoKind::Qua4;

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
        const NYM: usize = 4; // number of divisions in M
        const NYN: usize = 3; // number of divisions in N
        const NYL: usize = 2; // number of divisions for the extra height on the left side
        const NYG: usize = 4; // number of divisions for the gap
        const NYP: usize = 4; // number of divisions for the pile
        block1.set_ndiv(&[NYM, NYG]).unwrap();
        block2.set_ndiv(&[NYM, NYP]).unwrap();
        block3.set_ndiv(&[NYM, NYL]).unwrap();
        block4.set_ndiv(&[1, NYG]).unwrap();
        block5.set_ndiv(&[NYN, NYG]).unwrap();
        block6.set_ndiv(&[NYN, NYP]).unwrap();

        // generate mesh
        let mesh1 = block1.subdivide(kind).unwrap();
        let mesh2 = block2.subdivide(kind).unwrap();
        let mesh3 = block3.subdivide(kind).unwrap();
        let mesh4 = block4.subdivide(kind).unwrap();
        let mesh5 = block5.subdivide(kind).unwrap();
        let mesh6 = block6.subdivide(kind).unwrap();
        let mesh = join_meshes(&[&mesh1, &mesh2, &mesh3, &mesh4, &mesh5, &mesh6]).unwrap();
        mesh.check_all();

        // draw figure
        let mut fig = Figure::new();
        fig.show_point_ids(true)
            .show_cell_ids(true)
            .show_cell_att(false)
            .size(800.0, 800.0)
            .draw(&mesh, &format!("/tmp/pmsim/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
