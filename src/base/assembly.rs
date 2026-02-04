use crate::StrError;
use russell_lab::Matrix;
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};

/// Checks the symmetry of a local matrix
pub fn check_symmetry_of_local_matrix(kke: &Matrix, tol: f64) -> Result<(), StrError> {
    let n_equation_local = kke.dims().0;
    for l in 0..n_equation_local {
        for ll in (l + 1)..n_equation_local {
            if f64::abs(kke.get(l, ll) - kke.get(ll, l)) > tol {
                return Err("local matrix is not symmetric");
            }
        }
    }
    Ok(())
}

/// Assembles the local K matrix into its global counterpart for the Lagrange Multipliers Method (LMM)
pub fn assemble_matrix_lmm(kk: &mut CooMatrix, kke: &Matrix, local_to_global: &[usize]) -> Result<(), StrError> {
    let sym = kk.get_info().3;
    let n_equation_local = local_to_global.len();
    match sym {
        Sym::YesLower => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                for ll in 0..n_equation_local {
                    let gg = local_to_global[ll];
                    if g >= gg {
                        kk.put(g, gg, kke.get(l, ll)).unwrap();
                    }
                }
            }
        }
        Sym::YesUpper => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                for ll in 0..n_equation_local {
                    let gg = local_to_global[ll];
                    if g <= gg {
                        kk.put(g, gg, kke.get(l, ll)).unwrap();
                    }
                }
            }
        }
        Sym::YesFull | Sym::No => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                for ll in 0..n_equation_local {
                    let gg = local_to_global[ll];
                    kk.put(g, gg, kke.get(l, ll)).unwrap();
                }
            }
        }
    }
    Ok(())
}

/// Increments the number of non-zeros on the global K-bar and K-check matrices for the System Partitioning Strategy (SPS)
///
/// The Systems Partitioning Strategy (SPS) considers the following partitioning:
///
/// ```text
/// ┌       ┐ ┌   ┐   ┌   ┐
/// │ K̄   Ǩ │ │ ̄a │   │ f̄ │
/// │       │ │   │ = │   │
/// │ Ḵ   ̰K │ │ ǎ │   │ f̌ │
/// └       ┘ └   ┘   └   ┘
///     K       a       f
/// ```
pub fn add_nnz_sps(
    nnz_kk_bar: &mut usize,
    nnz_kk_check: &mut usize,
    sym: Sym,
    local_to_global: &[usize],
    eq_handler: &EquationHandler,
) {
    let n_equation_local = local_to_global.len();
    match sym {
        Sym::YesLower => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            if g >= gg {
                                *nnz_kk_bar += 1;
                            }
                        } else {
                            *nnz_kk_check += 1;
                        }
                    }
                }
            }
        }
        Sym::YesUpper => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            if g <= gg {
                                *nnz_kk_bar += 1;
                            }
                        } else {
                            *nnz_kk_check += 1;
                        }
                    }
                }
            }
        }
        Sym::YesFull | Sym::No => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            *nnz_kk_bar += 1;
                        } else {
                            *nnz_kk_check += 1;
                        }
                    }
                }
            }
        }
    }
}

/// Assembles the local K̄ and Ǩ matrices into their global counterparts for the System Partitioning Strategy (SPS)
///
/// The Systems Partitioning Strategy (SPS) considers the following partitioning:
///
/// ```text
/// ┌       ┐ ┌   ┐   ┌   ┐
/// │ K̄   Ǩ │ │ ̄a │   │ f̄ │
/// │       │ │   │ = │   │
/// │ Ḵ   ̰K │ │ ǎ │   │ f̌ │
/// └       ┘ └   ┘   └   ┘
///     K       a       f
/// ```
pub fn assemble_matrix_sps(
    kk_bar: &mut CooMatrix,
    kk_check: &mut CooMatrix,
    kke: &Matrix,
    local_to_global: &[usize],
    eq_handler: &EquationHandler,
) -> Result<(), StrError> {
    let sym = kk_bar.get_info().3;
    let n_equation_local = local_to_global.len();
    match sym {
        Sym::YesLower => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    let i = eq_handler.iu(g);
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            if g >= gg {
                                let j = eq_handler.iu(gg);
                                kk_bar.put(i, j, kke.get(l, ll)).unwrap();
                            }
                        } else {
                            let j = eq_handler.ip(gg);
                            kk_check.put(i, j, kke.get(l, ll)).unwrap();
                        }
                    }
                }
            }
        }
        Sym::YesUpper => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    let i = eq_handler.iu(g);
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            if g <= gg {
                                let j = eq_handler.iu(gg);
                                kk_bar.put(i, j, kke.get(l, ll)).unwrap();
                            }
                        } else {
                            let j = eq_handler.ip(gg);
                            kk_check.put(i, j, kke.get(l, ll)).unwrap();
                        }
                    }
                }
            }
        }
        Sym::YesFull | Sym::No => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if eq_handler.is_unknown(g) {
                    let i = eq_handler.iu(g);
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if eq_handler.is_unknown(gg) {
                            let j = eq_handler.iu(gg);
                            kk_bar.put(i, j, kke.get(l, ll)).unwrap();
                        } else {
                            let j = eq_handler.ip(gg);
                            kk_check.put(i, j, kke.get(l, ll)).unwrap();
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
#[cfg(test)]
mod tests {
    use super::assemble_matrix_lmm;
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};

    #[test]
    fn assemble_matrix_lmm_works_unsymmetric_1() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates marker
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let nnz = neq * neq;
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::No).unwrap();
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
        let mut ignore = vec![false; neq];
        ignore[2] = true;
        let tol = Some(1e-12);
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        #[rustfmt::skip]
        let correct = &[
            [10.0,     11.0, /*prescribed*/ 0.0,      0.0,   14.0], // 0
            [10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0], // 1
            [ 0.0,      0.0, /*prescribed*/ 0.0,      0.0,    0.0], // 2 (all prescribed)
            [ 0.0, 312100.0, /*prescribed*/ 0.0, 332300.0, 2400.0], // 3
            [10.0,   2111.0, /*prescribed*/ 0.0,   2300.0, 2414.0], // 4
        ];
        mat_approx_eq(&mat, correct, 1e-15);
    }

    #[test]
    fn assemble_matrix_lmm_works_symmetric_1() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates marker
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let nnz = neq * neq;
        const WRONG: f64 = 1234.0;
        #[rustfmt::skip]
        let k0_wrong = Matrix::from(&[
            [1.0,   4.0, 6.0],
            [4.0,   2.0, 5.0],
            [WRONG, 5.0, 3.0],
        ]);
        #[rustfmt::skip]
        let k0 = Matrix::from(&[
            [1.0, 4.0, 6.0],
            [4.0, 2.0, 5.0],
            [6.0, 5.0, 3.0],
        ]);
        #[rustfmt::skip]
        let k1 = Matrix::from(&[
            [100.0, 400.0, 600.0],
            [400.0, 200.0, 500.0],
            [600.0, 500.0, 300.0],
        ]);
        #[rustfmt::skip]
        let k2 = Matrix::from(&[
            [1000.0, 4000.0, 6000.0],
            [4000.0, 2000.0, 5000.0],
            [6000.0, 5000.0, 3000.0],
        ]);

        let mut ignore = vec![false; neq];
        ignore[2] = true;
        let tol = Some(1e-12);

        #[rustfmt::skip]
        let kk_correct = &[
            [1.0,      4.0, /*prescribed*/ 0.0,      0.0,    6.0], // 0
            [4.0,   1102.0, /*prescribed*/ 0.0,   6400.0,  605.0], // 1
            [0.0,      0.0, /*prescribed*/ 0.0,      0.0,    0.0], // 2 (all prescribed)
            [0.0,   6400.0, /*prescribed*/ 0.0,   3200.0,  500.0], // 3
            [6.0,    605.0, /*prescribed*/ 0.0,    500.0,  303.0], // 4
        ];

        // capture non-symmetric local matrix
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::YesLower).unwrap();
        assert_eq!(
            assemble_matrix(&mut kk, &k0_wrong, &l2g[0], &ignore, tol).err(),
            Some("local matrix is not symmetric")
        );

        // lower
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // upper
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::YesUpper).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // full
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::YesFull).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);
    }

    #[test]
    fn assemble_matrix_lmm_works_unsymmetric_2() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates marker
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let nnz_sup = 3 * 3 * 3;
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::No).unwrap();
        #[rustfmt::skip]
        let k0 = Matrix::from(&[
            [1.0, 4.0, 6.0],
            [7.0, 2.0, 5.0],
            [9.0, 8.0, 3.0],
        ]);
        #[rustfmt::skip]
        let k1 = Matrix::from(&[
            [100.0, 400.0, 600.0],
            [700.0, 200.0, 500.0],
            [900.0, 800.0, 300.0],
        ]);
        #[rustfmt::skip]
        let k2 = Matrix::from(&[
            [1000.0, 4000.0, 6000.0],
            [7000.0, 2000.0, 5000.0],
            [9000.0, 8000.0, 3000.0],
        ]);
        let ignore = vec![false; neq];
        let tol = Some(1e-12);
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        #[rustfmt::skip]
        let correct = &[
            [1.0,   4.0,     0.0,    0.0,   6.0], // 0
            [7.0, 1102.0, 4000.0, 6400.0, 605.0], // 1
            [0.0, 7000.0, 2000.0, 5000.0,   0.0], // 2
            [0.0, 9700.0, 8000.0, 3200.0, 500.0], // 3
            [9.0,  908.0,    0.0,  800.0, 303.0], // 4
        ];
        mat_approx_eq(&mat, correct, 1e-15);
    }

    #[test]
    fn assemble_matrix_lmm_works_symmetric_2() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates marker
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let nnz_sup = 3 * 3 * 3;
        const WRONG: f64 = 1234.0;
        #[rustfmt::skip]
        let k0_wrong = Matrix::from(&[
            [1.0,   4.0, 6.0],
            [4.0,   2.0, 5.0],
            [WRONG, 5.0, 3.0],
        ]);
        #[rustfmt::skip]
        let k0 = Matrix::from(&[
            [1.0, 4.0, 6.0],
            [4.0, 2.0, 5.0],
            [6.0, 5.0, 3.0],
        ]);
        #[rustfmt::skip]
        let k1 = Matrix::from(&[
            [100.0, 400.0, 600.0],
            [400.0, 200.0, 500.0],
            [600.0, 500.0, 300.0],
        ]);
        #[rustfmt::skip]
        let k2 = Matrix::from(&[
            [1000.0, 4000.0, 6000.0],
            [4000.0, 2000.0, 5000.0],
            [6000.0, 5000.0, 3000.0],
        ]);

        #[rustfmt::skip]
        let kk_correct = &[
            [1.0,    4.0,    0.0,    0.0,   6.0], // 0
            [4.0, 1102.0, 4000.0, 6400.0, 605.0], // 1
            [0.0, 4000.0, 2000.0, 5000.0,   0.0], // 2
            [0.0, 6400.0, 5000.0, 3200.0, 500.0], // 3
            [6.0,  605.0,    0.0,  500.0, 303.0], // 4
        ];

        let ignore = vec![false; neq];
        let tol = Some(1e-12);

        // capture non-symmetric local matrix
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::YesLower).unwrap();
        assert_eq!(
            assemble_matrix(&mut kk, &k0_wrong, &l2g[0], &ignore, tol).err(),
            Some("local matrix is not symmetric")
        );

        // lower
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // upper
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::YesUpper).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // full
        let mut kk = CooMatrix::new(neq, neq, nnz_sup, Sym::YesFull).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &ignore, tol).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &ignore, tol).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);
    }
}
*/
