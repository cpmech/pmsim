use super::{ElementDofsMap, AllDofs};
use crate::StrError;
use gemlab::mesh::Cell;
use russell_lab::{Matrix, Vector};
use russell_sparse::{CooMatrix, Sym};

/// Computes local-to-global maps needed for the assembly process
pub fn compute_local_to_global(emap: &ElementDofsMap, eqs: &AllDofs, cell: &Cell) -> Result<Vec<usize>, StrError> {
    let info = emap.get(cell)?;
    let mut local_to_global = vec![0; info.n_equation];
    for m in 0..cell.points.len() {
        for (dof, local) in &info.dofs[m] {
            let global = eqs.eq(cell.points[m], *dof)?;
            local_to_global[*local] = global;
        }
    }
    Ok(local_to_global)
}

/// Assembles local vector into global vector
///
/// # Output
///
/// * `rr` -- is the global vector R with length = `n_equation`
///
/// # Input
///
/// * `phi` -- is the local vector with length = `n_equation_local`
/// * `cell_id` -- is the ID of the cell adding the contribution to R
/// * `local_to_global` -- is an array holding all equation numbers.
/// * `ignore` -- (n_equation) holds the equation numbers to be ignored in the assembly process;
///   i.e., it allows the generation of the reduced system. For example, the equations corresponding
///   to the prescribed DOFs must not be assembled in the reduced system.
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
pub fn assemble_vector(rr: &mut Vector, phi: &Vector, local_to_global: &[usize], ignore: &[bool]) {
    let n_equation_local = phi.dim();
    for l in 0..n_equation_local {
        let g = local_to_global[l];
        if !ignore[g] {
            rr[g] += phi[l];
        }
    }
}

/// Assembles local matrix into global matrix
///
/// # Output
///
/// * `kk` -- is the global square matrix K with dims = (`n_equation`,`n_equation`)
///
/// # Input
///
/// * `kke` -- is the local square matrix Ke with dims = (`n_equation_local`,`n_equation_local`)
/// * `cell_id` -- is the ID of the cell adding the contribution to K
/// * `local_to_global` -- is an nested holding all equation numbers.
/// * `ignore` -- (n_equation) holds the equation numbers to be ignored in the assembly process;
///   i.e., it allows the generation of the reduced system. For example, the equations corresponding
///   to the prescribed DOFs must not be assembled in the reduced system.
///
/// # Note
///
/// Will check the symmetry of the local matrix if the global matrix has the symmetry flag enabled.
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
pub fn assemble_matrix(
    kk: &mut CooMatrix,
    kke: &Matrix,
    local_to_global: &[usize],
    ignore: &[bool],
    symmetry_check_tolerance: Option<f64>,
) -> Result<(), StrError> {
    let n_equation_local = kke.dims().0;
    // check symmetry of local matrices
    let sym = kk.get_info().3;
    let symmetric = sym != Sym::No;
    if symmetric {
        if let Some(tol) = symmetry_check_tolerance {
            for l in 0..n_equation_local {
                for ll in (l + 1)..n_equation_local {
                    if f64::abs(kke.get(l, ll) - kke.get(ll, l)) > tol {
                        return Err("local matrix is not symmetric");
                    }
                }
            }
        }
    }
    // assemble
    match sym {
        Sym::YesLower => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !ignore[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !ignore[gg] {
                            if g >= gg {
                                kk.put(g, gg, kke.get(l, ll)).unwrap();
                            }
                        }
                    }
                }
            }
        }
        Sym::YesUpper => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !ignore[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !ignore[gg] {
                            if g <= gg {
                                kk.put(g, gg, kke.get(l, ll)).unwrap();
                            }
                        }
                    }
                }
            }
        }
        Sym::YesFull | Sym::No => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !ignore[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !ignore[gg] {
                            kk.put(g, gg, kke.get(l, ll)).unwrap();
                        }
                    }
                }
            }
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{assemble_matrix, assemble_vector};
    use crate::base::{compute_local_to_global, Attributes, Elem, ElementDofsMap, AllDofs};
    use crate::base::{ParamBeam, ParamPorousLiq, ParamPorousSldLiq, ParamSolid};
    use gemlab::{mesh::Samples, shapes::GeoKind};
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};

    #[test]
    fn compute_local_to_global_handles_errors() {
        let mut mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let amap = Attributes::from([(1, Elem::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        let eqs = AllDofs::new(&mesh, &emap).unwrap();
        mesh.cells[0].kind = GeoKind::Qua4; // never do this!
        assert_eq!(
            compute_local_to_global(&emap, &eqs, &mesh.cells[0]).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
        mesh.cells[0].kind = GeoKind::Tri3;
        mesh.cells[0].points[0] = 100; // never do this!
        assert_eq!(
            compute_local_to_global(&emap, &eqs, &mesh.cells[0]).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
    }

    #[test]
    fn compute_local_to_global_works() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let amap = Attributes::from([(1, Elem::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        let eqs = AllDofs::new(&mesh, &emap).unwrap();
        let l2g0 = compute_local_to_global(&emap, &eqs, &mesh.cells[0]).unwrap();
        let l2g1 = compute_local_to_global(&emap, &eqs, &mesh.cells[1]).unwrap();
        let l2g2 = compute_local_to_global(&emap, &eqs, &mesh.cells[2]).unwrap();
        assert_eq!(l2g0, &[0, 1, 2, 3, 8, 9]);
        assert_eq!(l2g1, &[2, 3, 6, 7, 8, 9]);
        assert_eq!(l2g2, &[2, 3, 4, 5, 6, 7]);

        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let amap = Attributes::from([(1, Elem::PorousLiq(p)), (2, Elem::PorousLiq(p))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        let eqs = AllDofs::new(&mesh, &emap).unwrap();
        let l2g0 = compute_local_to_global(&emap, &eqs, &mesh.cells[0]).unwrap();
        let l2g1 = compute_local_to_global(&emap, &eqs, &mesh.cells[1]).unwrap();
        let l2g2 = compute_local_to_global(&emap, &eqs, &mesh.cells[2]).unwrap();
        assert_eq!(l2g0, &[0, 1, 3]);
        assert_eq!(l2g1, &[2, 3, 1]);
        assert_eq!(l2g2, &[1, 4, 5, 2]);

        // 8------7------6._
        // |       [3](3)|  '-.5
        // |  [0]        |     '-._
        // 9  (1)       10  [1]    '4
        // |             |  (2)  .-'
        // |       [2](3)|   _.3'
        // 0------1------2.-'
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let amap = Attributes::from([(1, Elem::PorousSldLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Beam(p3))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        let eqs = AllDofs::new(&mesh, &emap).unwrap();
        let l2g0 = compute_local_to_global(&emap, &eqs, &mesh.cells[0]).unwrap();
        let l2g1 = compute_local_to_global(&emap, &eqs, &mesh.cells[1]).unwrap();
        let l2g2 = compute_local_to_global(&emap, &eqs, &mesh.cells[2]).unwrap();
        let l2g3 = compute_local_to_global(&emap, &eqs, &mesh.cells[3]).unwrap();
        assert_eq!(
            l2g0,
            &[0, 1, 5, 6, 15, 16, 21, 22, 3, 4, 26, 27, 19, 20, 24, 25, 2, 8, 18, 23]
        );
        assert_eq!(l2g1, &[5, 6, 11, 12, 15, 16, 9, 10, 13, 14, 26, 27]);
        assert_eq!(l2g2, &[5, 6, 7, 26, 27, 28]);
        assert_eq!(l2g3, &[26, 27, 28, 15, 16, 17]);
    }

    #[test]
    fn assemble_vector_works_1() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
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
        let mut ignore = vec![false; neq];
        ignore[2] = true;
        assemble_vector(&mut ff, &f0, &l2g[0], &ignore);
        assemble_vector(&mut ff, &f1, &l2g[1], &ignore);
        assemble_vector(&mut ff, &f2, &l2g[2], &ignore);
        assert_eq!(ff.as_data(), &[10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0]);
    }

    #[test]
    fn assemble_vector_works_2() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let l2g = vec![vec![0, 1, 4], vec![1, 3, 4], vec![1, 2, 3]];
        let neq = 5;
        let mut ff = Vector::new(neq);
        let f0 = Vector::from(&[1.0, 2.0, 3.0]);
        let f1 = Vector::from(&[10.0, 20.0, 30.0]);
        let f2 = Vector::from(&[100.0, 200.0, 300.0]);
        let ignore = vec![false; neq];
        assemble_vector(&mut ff, &f0, &l2g[0], &ignore);
        assemble_vector(&mut ff, &f1, &l2g[1], &ignore);
        assemble_vector(&mut ff, &f2, &l2g[2], &ignore);
        assert_eq!(ff.as_data(), &[1.0, 112.0, 200.0, 320.0, 33.0]);
    }

    #[test]
    fn assemble_matrix_works_unsymmetric_1() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
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
    fn assemble_matrix_works_symmetric_1() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
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
    fn assemble_matrix_works_unsymmetric_2() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
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
    fn assemble_matrix_works_symmetric_2() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
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
