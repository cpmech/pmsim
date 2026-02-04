use crate::StrError;
use russell_lab::Matrix;
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};

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
