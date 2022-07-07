use super::{Equation, NodalDofs};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::Vector;

/// Implements functions to perform the assembly process
pub struct Assemble {
    /// Holds all local-to-global mappings (ncell x n_local_equation)
    pub local_to_global: Vec<Vec<usize>>,

    /// Records what equations are prescribed (nequation)
    pub prescribed: Vec<bool>,
}

impl Assemble {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, nodal_dofs: &NodalDofs, equation: &Equation) -> Result<Self, StrError> {
        let mut local_to_global = vec![Vec::new(); mesh.cells.len()];
        for cell in &mesh.cells {
            let dofs = match nodal_dofs.get(cell.attribute_id, cell.kind) {
                Some(d) => d,
                None => return Err("cannot find nodal DOFs for given (CellAttributeId,GeoKind)"),
            };
            for m in 0..cell.points.len() {
                for dof in &dofs[m] {
                    let eid = match equation.id(cell.points[m], *dof) {
                        Some(e) => e,
                        None => return Err("cannot find equation id for given (PointId,DOF)"),
                    };
                    local_to_global[cell.id].push(eid);
                }
            }
        }
        Ok(Assemble {
            local_to_global,
            prescribed: vec![false; equation.len()],
        })
    }

    /// Assembles local vector into global vector
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
    use crate::base::{AttrElement, Element, Equation, NodalDofs, SampleMeshes};
    use crate::StrError;
    use russell_lab::Vector;

    #[test]
    fn assemble_works() -> Result<(), StrError> {
        //       {2}4---.__
        //         / \     `--.___3{3}   [#] indicates id
        //        /   \          / \     (#) indicates attribute_id
        //       /     \  [1]   /   \    {#} indicates equation id
        //      /  [0]  \ (1)  / [2] \
        //     /   (1)   \    /  (1)  \
        // {0}0---.__     \  /      ___2{4}
        //           `--.__\/__.---'
        //               {1}1
        let mesh = SampleMeshes::three_tri3();
        let attr_element = AttrElement::from([(1, Element::PorousLiq)]);
        let nodal_dofs = NodalDofs::new(&mesh, &attr_element).unwrap();
        let equation = Equation::new(&mesh, &nodal_dofs).unwrap();
        let assemble = Assemble::new(&mesh, &nodal_dofs, &equation)?;
        assert_eq!(assemble.local_to_global, &[[0, 1, 2], [1, 3, 2], [1, 4, 3]]);

        let mut global = Vector::new(equation.len());
        let local_0 = Vector::from(&[/*    */ 10.0, /*    */ 11.0, /*    */ 14.0]);
        let local_1 = Vector::from(&[/*  */ 2100.0, /*  */ 2300.0, /*  */ 2400.0]);
        let local_2 = Vector::from(&[/**/ 310000.0, /**/ 320000.0, /**/ 330000.0]);
        assemble.vector(&mut global, &local_0, 0);
        assemble.vector(&mut global, &local_1, 1);
        assemble.vector(&mut global, &local_2, 2);
        assert_eq!(global.as_data(), &[10.0, 312111.0, 2414.0, 332300.0, 320000.0]);
        Ok(())
    }
}
