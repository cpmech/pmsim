use super::{ElementDofsMap, Equations};
use crate::StrError;
use gemlab::mesh::Cell;
use russell_lab::{Matrix, Vector};
use russell_sparse::{CooMatrix, Sym};

const SYMMETRY_CHECK_TOLERANCE: f64 = 1e-12;

/// Computes local-to-global maps needed for the assembly process
pub fn compute_local_to_global(info: &ElementDofsMap, eqs: &Equations, cell: &Cell) -> Result<Vec<usize>, StrError> {
    let info = info.get(cell)?;
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
/// # Note
///
/// Will check the symmetry of the local matrix if the global matrix has the symmetry flag enabled.
///
/// # Panics
///
/// This function will panic if the indices are out-of-bounds
pub fn assemble_matrix(
    kk_global: &mut CooMatrix,
    kk_local: &Matrix,
    local_to_global: &[usize],
    prescribed: &[bool],
) -> Result<(), StrError> {
    let n_equation_local = kk_local.dims().0;
    // check symmetry of local matrices
    let sym = kk_global.get_info().3;
    let symmetric = sym != Sym::No;
    if symmetric {
        for l in 0..n_equation_local {
            for ll in (l + 1)..n_equation_local {
                if f64::abs(kk_local.get(l, ll) - kk_local.get(ll, l)) > SYMMETRY_CHECK_TOLERANCE {
                    return Err("local matrix is not symmetric");
                }
            }
        }
    }
    // assemble
    match sym {
        Sym::YesLower => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !prescribed[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !prescribed[gg] {
                            if g >= gg {
                                kk_global.put(g, gg, kk_local.get(l, ll)).unwrap();
                            }
                        }
                    }
                }
            }
        }
        Sym::YesUpper => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !prescribed[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !prescribed[gg] {
                            if g <= gg {
                                kk_global.put(g, gg, kk_local.get(l, ll)).unwrap();
                            }
                        }
                    }
                }
            }
        }
        Sym::YesFull | Sym::No => {
            for l in 0..n_equation_local {
                let g = local_to_global[l];
                if !prescribed[g] {
                    for ll in 0..n_equation_local {
                        let gg = local_to_global[ll];
                        if !prescribed[gg] {
                            kk_global.put(g, gg, kk_local.get(l, ll)).unwrap();
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
    use crate::base::{compute_local_to_global, Attributes, Element, ElementDofsMap, Equations};
    use crate::base::{ParamBeam, ParamPorousLiq, ParamPorousSldLiq, ParamSolid};
    use gemlab::{mesh::Samples, shapes::GeoKind};
    use russell_lab::{mat_approx_eq, Matrix, Vector};
    use russell_sparse::{CooMatrix, Sym};

    #[test]
    fn compute_local_to_global_handles_errors() {
        let mut mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
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
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
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
        let att = Attributes::from([(1, Element::PorousLiq(p)), (2, Element::PorousLiq(p))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
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
        let att = Attributes::from([
            (1, Element::PorousSldLiq(p1)),
            (2, Element::Solid(p2)),
            (3, Element::Beam(p3)),
        ]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
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
    fn assemble_vector_works() {
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
        let mut prescribed = vec![false; neq];
        prescribed[2] = true;
        assemble_vector(&mut ff, &f0, &l2g[0], &prescribed);
        assemble_vector(&mut ff, &f1, &l2g[1], &prescribed);
        assemble_vector(&mut ff, &f2, &l2g[2], &prescribed);
        assert_eq!(ff.as_data(), &[10.0, 312111.0, /*prescribed*/ 0.0, 332300.0, 2414.0]);
    }

    #[test]
    fn assemble_matrix_works_unsymmetric() {
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
        let mut prescribed = vec![false; neq];
        prescribed[2] = true;
        assemble_matrix(&mut kk, &k0, &l2g[0], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &prescribed).unwrap();
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
    fn assemble_matrix_works_symmetric() {
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

        let mut prescribed = vec![false; neq];
        prescribed[2] = true;

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
            assemble_matrix(&mut kk, &k0_wrong, &l2g[0], &prescribed).err(),
            Some("local matrix is not symmetric")
        );

        // lower
        assemble_matrix(&mut kk, &k0, &l2g[0], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &prescribed).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // upper
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::YesUpper).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &prescribed).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);

        // full
        let mut kk = CooMatrix::new(neq, neq, nnz, Sym::YesFull).unwrap();
        assemble_matrix(&mut kk, &k0, &l2g[0], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k1, &l2g[1], &prescribed).unwrap();
        assemble_matrix(&mut kk, &k2, &l2g[2], &prescribed).unwrap();
        let mat = kk.as_dense();
        mat_approx_eq(&mat, kk_correct, 1e-15);
    }
}
