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
        const P: f64 = 6.0; // pile wall height
        const R: f64 = P + 2.25; // right height
        const L: f64 = R + 2.0; // left height
        const M: f64 = 1.0 * L; // left distance (it should be > 1.5 L)
        const G: f64 = 0.5; // gap (thickness of sheet pile wall)
        const N: f64 = 5.5 / 2.0; // right distance
        const NX1: usize = 4; // number of divisions along x
        const NX2: usize = 2; // number of divisions along x
        const NY1: usize = 3; // number of divisions along y
        const NY2: usize = 2; // number of divisions along y
        const NY3: usize = 1; // number of divisions along y
        let kind = GeoKind::Qua4;

        // generate mesh
        let mut block1 = Block::new(&[[0.0, 0.0], [M, 0.0], [M, P], [0.0, P]]).unwrap();
        let mut block2 = Block::new(&[[0.0, P], [M, P], [M, L], [0.0, L]]).unwrap();
        let mut block3 = Block::new(&[[M, 0.0], [M + G, 0.0], [M + G, P], [M, P]]).unwrap();
        let mut block4 = Block::new(&[[M + G, 0.0], [M + G + N, 0.0], [M + G + N, P], [M + G, P]]).unwrap();
        let mut block5 = Block::new(&[[M + G, P], [M + G + N, P], [M + G + N, R], [M + G, R]]).unwrap();

        block1.set_ndiv(&[NX1, NY1]).unwrap();
        block2.set_ndiv(&[NX1, NY2]).unwrap();
        block3.set_ndiv(&[1, NY1]).unwrap();
        block4.set_ndiv(&[NX2, NY1]).unwrap();
        block5.set_ndiv(&[NX2, NY3]).unwrap();

        let mesh1 = block1.subdivide(kind).unwrap();
        let mesh2 = block2.subdivide(kind).unwrap();
        let mesh3 = block3.subdivide(kind).unwrap();
        let mesh4 = block4.subdivide(kind).unwrap();
        let mesh5 = block5.subdivide(kind).unwrap();
        let mesh = join_meshes(&[&mesh1, &mesh2, &mesh3, &mesh4, &mesh5]).unwrap();
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
