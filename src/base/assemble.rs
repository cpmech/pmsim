use super::LocalToGlobal;
use gemlab::mesh::CellId;
use russell_lab::{Matrix, Vector};
use russell_sparse::SparseTriplet;

/// Assembles local vector into global vector
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
#[inline]
pub fn assemble_vector(ff: &mut Vector, f: &Vector, cell_id: CellId, l2g: &LocalToGlobal, prescribed: &Vec<bool>) {
    for l in 0..f.dim() {
        let g = l2g[cell_id][l];
        if !prescribed[g] {
            ff[g] += f[l];
        }
    }
}

/// Assembles local matrix into global matrix
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
#[inline]
pub fn assemble_matrix(
    kk: &mut SparseTriplet,
    k: &Matrix,
    cell_id: CellId,
    l2g: &LocalToGlobal,
    prescribed: &Vec<bool>,
) {
    let nlocal = k.dims().0;
    for l in 0..nlocal {
        let g = l2g[cell_id][l];
        if !prescribed[g] {
            for ll in 0..nlocal {
                let gg = l2g[cell_id][ll];
                if !prescribed[gg] {
                    kk.put(g, gg, k[l][ll]).unwrap();
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{assemble_matrix, assemble_vector};
    use russell_lab::{Matrix, Vector};
    use russell_sparse::{SparseTriplet, Symmetry};

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
        assemble_vector(&mut ff, &f0, 0, &l2g, &prescribed);
        assemble_vector(&mut ff, &f1, 1, &l2g, &prescribed);
        assemble_vector(&mut ff, &f2, 2, &l2g, &prescribed);
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
        assemble_matrix(&mut kk, &k0, 0, &l2g, &prescribed);
        assemble_matrix(&mut kk, &k1, 1, &l2g, &prescribed);
        assemble_matrix(&mut kk, &k2, 2, &l2g, &prescribed);
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
