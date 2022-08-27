use super::{BcsEssential, DofNumbers};
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Generates the array indicating the prescribed DOFs
///
/// Returns the array `prescribed` indicating which DOFs (equations) are prescribed.
/// The length of the array is equal to the total number of DOFs which is
/// equal to the total number of equations `n_equation`.
pub fn gen_prescribed_array(dn: &DofNumbers, ebc: &BcsEssential) -> Result<Vec<bool>, StrError> {
    let mut prescribed = vec![false; dn.n_equation];
    for (point_id, dof) in ebc.all.keys() {
        match dn.point_dofs[*point_id].get(dof) {
            Some(eq) => prescribed[*eq] = true,
            None => return Err("EBC dof is not present in point_dofs array"),
        }
    }
    Ok(prescribed)
}

/// Assembles local vector into global vector
///
/// # Output
///
/// * `rr_global` -- is the global vector R with length = `n_equation`
///
/// # Input
///
/// * `r_local` -- is the local vector r with length = `n_equation_local`
/// * `cell_id` -- is the ID of the cell adding the contribution to R
/// * `local_to_global` -- is an array holding all equation numbers.
/// * `prescribed` -- tells whether a global equation number has prescribed
///   DOF or not. Its length is equal to the total number of DOFs `n_equation`.
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
#[inline]
pub fn assemble_vector(rr_global: &mut Vector, r_local: &Vector, local_to_global: &[usize], prescribed: &[bool]) {
    let n_equation_local = r_local.dim();
    for l in 0..n_equation_local {
        let g = local_to_global[l];
        if !prescribed[g] {
            rr_global[g] += r_local[l];
        }
    }
}

/// Assembles local matrix into global matrix
///
/// # Output
///
/// * `kk_global` -- is the global square matrix K with dims = (`n_equation`,`n_equation`)
///
/// # Input
///
/// * `kk_local` -- is the local square matrix K with dims = (`n_equation_local`,`n_equation_local`)
/// * `cell_id` -- is the ID of the cell adding the contribution to K
/// * `local_to_global` -- is an nested holding all equation numbers.
/// * `prescribed` -- tells whether a global equation number has prescribed
///   DOF or not. Its length is equal to the total number of DOFs `n_equation`.
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
#[inline]
pub fn assemble_matrix(
    kk_global: &mut SparseTriplet,
    kk_local: &Matrix,
    local_to_global: &[usize],
    prescribed: &[bool],
) {
    let n_equation_local = kk_local.dims().0;
    for l in 0..n_equation_local {
        let g = local_to_global[l];
        if !prescribed[g] {
            for ll in 0..n_equation_local {
                let gg = local_to_global[ll];
                if !prescribed[gg] {
                    kk_global.put(g, gg, kk_local[l][ll]).unwrap();
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{assemble_matrix, assemble_vector, gen_prescribed_array};
    use crate::base::{BcsEssential, Dof, DofNumbers, Element, SampleParams};
    use gemlab::mesh::Samples;
    use russell_lab::{Matrix, Vector};
    use russell_sparse::{SparseTriplet, Symmetry};
    use std::collections::HashMap;

    #[test]
    fn gen_prescribed_array_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute_id
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_porous_liq();
        let elements = HashMap::from([(1, Element::PorousLiq(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let mut ebc = BcsEssential::new();
        let zero = |_| 0.0;
        ebc.at(&[0, 4], &[Dof::Pl], zero);
        let prescribed = gen_prescribed_array(&dn, &ebc).unwrap();
        assert_eq!(prescribed, &[true, false, false, false, true]);
        assert_eq!(zero(1.0), 0.0);

        ebc.at(&[3], &[Dof::T], zero);
        assert_eq!(
            gen_prescribed_array(&dn, &ebc).err(),
            Some("EBC dof is not present in point_dofs array")
        );

        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute_id
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let mut ebc = BcsEssential::new();
        ebc.at(&[0], &[Dof::Ux, Dof::Uy], zero);
        ebc.at(&[1, 2], &[Dof::Uy], zero);
        let prescribed = gen_prescribed_array(&dn, &ebc).unwrap();
        assert_eq!(
            prescribed,
            //   0     1      2     3      4     5      6      7      8      9
            &[true, true, false, true, false, true, false, false, false, false]
        );
    }

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
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let mut ff = Vector::new(neq);
        let f0 = Vector::from(&[/*    */ 10.0, /*    */ 11.0, /*    */ 14.0]);
        let f1 = Vector::from(&[/*  */ 2100.0, /*  */ 2300.0, /*  */ 2400.0]);
        let f2 = Vector::from(&[/**/ 310000.0, /**/ 320000.0, /**/ 330000.0]);
        let mut prescribed = vec![false; neq];
        prescribed[2] = true;
        assemble_vector(&mut ff, &f0, &l2g[0], &prescribed);
        assemble_vector(&mut ff, &f1, &l2g[1], &prescribed);
        assemble_vector(&mut ff, &f2, &l2g[2], &prescribed);
        assert_eq!(ff.as_data(), &[10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0]);
    }

    #[test]
    fn assemble_matrix_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute_id
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let mut kk = SparseTriplet::new(neq, neq, neq * neq, Symmetry::No).unwrap();
        #[rustfmt::skip]
        let k0 = Matrix::from(&[
            [10.0, 11.0, 14.0],
            [10.0, 11.0, 14.0],
            [10.0, 11.0, 14.0],
        ]);
        #[rustfmt::skip]
        let k1 = Matrix::from(&[
            [2100.0, 2300.0, 2400.0],
            [2100.0, 2300.0, 2400.0],
            [2100.0, 2300.0, 2400.0],
        ]);
        #[rustfmt::skip]
        let k2 = Matrix::from(&[
            [310000.0, 320000.0, 330000.0],
            [310000.0, 320000.0, 330000.0],
            [310000.0, 320000.0, 330000.0],
        ]);
        let mut prescribed = vec![false; neq];
        prescribed[2] = true;
        assemble_matrix(&mut kk, &k0, &l2g[0], &prescribed);
        assemble_matrix(&mut kk, &k1, &l2g[1], &prescribed);
        assemble_matrix(&mut kk, &k2, &l2g[2], &prescribed);
        let mut kk_mat = Matrix::new(5, 5);
        kk.to_matrix(&mut kk_mat).unwrap();
        #[rustfmt::skip]
        let correct = &[
            10.0,     11.0, /*prescribed*/ 0.0,      0.0,   14.0, // 0
            10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0, // 1
             0.0,      0.0, /*prescribed*/ 0.0,      0.0,    0.0, // 2 (all prescribed)
             0.0, 312100.0, /*prescribed*/ 0.0, 332300.0, 2400.0, // 3
            10.0,   2111.0, /*prescribed*/ 0.0,   2300.0, 2414.0, // 4
        ];
        assert_eq!(kk_mat.as_data(), correct);
    }
}
