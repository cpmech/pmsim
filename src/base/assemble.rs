#![allow(unused)]

use super::{Equation, LocalToGlobal, NodalDofs};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::Vector;

/// Assembles local vector into global vector
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
pub fn assemble_vector(
    global: &mut Vector,
    local: &Vector,
    cell_id: CellId,
    local_to_global: &LocalToGlobal,
    prescribed: &Vec<bool>,
) {
    let l2g = &local_to_global[cell_id];
    let nglobal = global.dim();
    let nlocal = local.dim();
    for l in 0..nlocal {
        let g = local_to_global[cell_id][l];
        if !prescribed[g] {
            global[g] += local[l];
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::assemble_vector;
    use crate::base::{AttrElement, Element, Equation, NodalDofs, SampleMeshes};
    use crate::StrError;
    use russell_lab::Vector;

    #[test]
    fn assemble_vector_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute_id
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let local_to_global = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let mut global = Vector::new(5);
        let local_0 = Vector::from(&[/*    */ 10.0, /*    */ 11.0, /*    */ 14.0]);
        let local_1 = Vector::from(&[/*  */ 2100.0, /*  */ 2300.0, /*  */ 2400.0]);
        let local_2 = Vector::from(&[/**/ 310000.0, /**/ 320000.0, /**/ 330000.0]);
        let mut prescribed = vec![false; global.dim()];
        prescribed[2] = true;
        assemble_vector(&mut global, &local_0, 0, &local_to_global, &prescribed);
        assemble_vector(&mut global, &local_1, 1, &local_to_global, &prescribed);
        assemble_vector(&mut global, &local_2, 2, &local_to_global, &prescribed);
        assert_eq!(
            global.as_data(),
            &[10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0]
        );
    }
}
