use super::{Dof, Equation};
use gemlab::mesh::{CellId, Mesh};
use russell_lab::Vector;

pub struct Assemble {
    /// Holds all local-to-global mappings (ncell x n_local_equation)
    pub local_to_global: Vec<Vec<usize>>,

    /// Records what equations are prescribed (nequation)
    pub prescribed: Vec<bool>,
}

impl Assemble {
    // cell_dofs [ncell][nnode]
    pub fn new(mesh: &Mesh, equation: &Equation, cell_dofs: &Vec<Vec<Dof>>) -> Self {
        let mut local_to_global = vec![Vec::new(); mesh.cells.len()];
        for cell in &mesh.cells {
            for m in 0..cell.points.len() {
                let eid = equation.id(cell.points[m], cell_dofs[cell.id][m]).unwrap();
                local_to_global[cell.id].push(eid);
            }
        }
        Assemble {
            local_to_global,
            prescribed: vec![false; equation.len()],
        }
    }

    pub fn vector(&self, global: &mut Vector, local: &Vector, cell_id: CellId) {
        for i in 0..local.dim() {
            let g = self.local_to_global[cell_id][i];
            if !self.prescribed[g] {
                global[g] += local[i];
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Assemble;
    use crate::base::{Dof, Equation, SampleMeshes};
    use russell_lab::Vector;

    #[test]
    fn new_and_vector_works() {
        //       {2}4---.__
        //         / \     `--.___3{3}
        //        /   \          / \
        //       /     \  [1]   /   \    {#} means equation id
        //      /  [0]  \      /     \
        //     /         \    /  [2]  \
        // {0}0---.__     \  /      ___2{4}
        //           `--.__\/__.---'
        //               {1}1
        let mesh = SampleMeshes::three_tri3();
        let mut equation = Equation::new(mesh.points.len());
        let cell_dofs = vec![
            vec![Dof::T, Dof::T, Dof::T],
            vec![Dof::T, Dof::T, Dof::T],
            vec![Dof::T, Dof::T, Dof::T],
        ];
        mesh.cells.iter().zip(&cell_dofs).for_each(|(cell, dofs)| {
            cell.points.iter().zip(dofs).for_each(|(point_id, dof)| {
                equation.activate(*point_id, *dof);
            });
        });
        let assemble = Assemble::new(&mesh, &equation, &cell_dofs);
        assert_eq!(assemble.local_to_global, &[[0, 1, 2], [1, 3, 2], [1, 4, 3]]);

        let mut global = Vector::new(equation.len());
        let local_0 = Vector::from(&[/*    */ 10.0, /*    */ 11.0, /*    */ 14.0]);
        let local_1 = Vector::from(&[/*  */ 2100.0, /*  */ 2300.0, /*  */ 2400.0]);
        let local_2 = Vector::from(&[/**/ 310000.0, /**/ 320000.0, /**/ 330000.0]);
        assemble.vector(&mut global, &local_0, 0);
        assemble.vector(&mut global, &local_1, 1);
        assemble.vector(&mut global, &local_2, 2);
        assert_eq!(global.as_data(), &[10.0, 312111.0, 2414.0, 332300.0, 320000.0]);
    }
}
