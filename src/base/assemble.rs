#![allow(unused)]

use super::Equation;
use gemlab::mesh::CellId;
use russell_lab::Vector;
use std::collections::HashSet;

pub struct Assemble {
    /// Holds all local-to-global mappings (ncell x n_local_equation)
    pub local_to_global: Vec<Vec<usize>>,
}

impl Assemble {
    pub fn new(ncell: usize) -> Self {
        Assemble {
            local_to_global: vec![Vec::new(); ncell],
        }
    }

    pub fn vector(&self, global: &mut Vector, local: &Vector, cell_id: CellId, prescribed: &Vec<bool>) {
        for i in 0..local.dim() {
            let g = self.local_to_global[cell_id][i];
            if !prescribed[g] {
                global[g] += local[i];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use russell_lab::Vector;

    use super::Assemble;
    use crate::base::SampleMeshes;

    #[test]
    fn new_works() {
        let assemble = Assemble::new(2);
        assert_eq!(assemble.local_to_global.len(), 2);
    }

    #[test]
    fn vector_works() {
        //         4---.__
        //        / \     `--.___3
        //       /   \          / \
        //      /     \  [1]   /   \
        //     /  [0]  \      /     \
        //    /         \    /  [2]  \
        //   0---.__     \  /      ___2
        //          `--.__\/__.---'
        //                 1
        let mesh = SampleMeshes::three_tri3();
        let mut assemble = Assemble::new(mesh.cells.len());
        for cell in &mesh.cells {
            for point_id in &cell.points {
                assemble.local_to_global[cell.id].push(*point_id);
            }
        }
        assert_eq!(assemble.local_to_global, &[[0, 1, 4], [1, 3, 4], [1, 2, 3]]);

        let mut global = Vector::new(5);
        let prescribed = vec![false; global.dim()];
        let local_0 = Vector::from(&[/*    */ 10.0, /*    */ 11.0, /*    */ 14.0]);
        let local_1 = Vector::from(&[/*  */ 2100.0, /*  */ 2300.0, /*  */ 2400.0]);
        let local_2 = Vector::from(&[/**/ 310000.0, /**/ 320000.0, /**/ 330000.0]);
        assemble.vector(&mut global, &local_0, 0, &prescribed);
        assemble.vector(&mut global, &local_1, 1, &prescribed);
        assemble.vector(&mut global, &local_2, 2, &prescribed);
        assert_eq!(global.as_data(), &[10.0, 312111.0, 320000.0, 332300.0, 2414.0]);
    }
}
