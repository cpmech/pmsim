use super::{ReferenceData, ReferenceDataType};
use crate::base::{Config, Dof, Schema, NZ_VON_MISES};
use crate::fem::{FemState, PostProc};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_tensor::SQRT_2;

/// Returns T or F for a boolean variable
fn b2s(flag: bool) -> String {
    if flag {
        "T".to_string()
    } else {
        "F".to_string()
    }
}

/// Queries whether A failed to compare with B or not
///
/// Returns `(fail, diff)`
fn query_failed(a: f64, b: f64, tol: f64, verbose: usize) -> (bool, f64) {
    let diff = f64::abs(a - b);
    let fail = diff > tol;
    if verbose == 1 {
        let mrk = if fail { "❌" } else { "➖" };
        print!("{:15.6e}{} ", diff, mrk);
    } else if verbose == 2 {
        let mrk = if fail { "❌" } else { "➖" };
        print!("{:9.2e} vs {:9.2e}({:9.2e}{}) ", a, b, diff, mrk);
    }
    (fail, diff)
}

/// Queries whether A equals B
///
/// Returns `fail`
fn query_failed_bool(a: bool, b: bool, verbose: usize) -> bool {
    let fail = a != b;
    if verbose > 0 {
        let mrk = if fail { "❌" } else { "➖" };
        print!("{} vs {} {}    ", b2s(a), b2s(b), mrk);
    }
    fail
}

/// Compares the FEM results (displacement, stress, strain) against reference data
///
/// # Input
///
/// * `mesh` -- The mesh
/// * `res_path` -- The full path to the results files; e.g., "/tmp/pmsim/simulation.json"
/// * `ref_type` -- The type (origin) of the reference data
/// * `ref_path` -- The full path of the file with the reference results
/// * `tol_displacement` -- A tolerance to compare displacements
/// * `tol_stress` -- A tolerance to compare stresses
/// * `verbose` -- Enables the verbose mode:
///   - 0 => no output
///   - 1 => shows error
///   - 2 => shows values and error
/// * `eps_bar_p` -- Tells this function to check the accumulated plastic strain (eps_bar_p).
///   The values in the tuple are `(index_in_z_set, conversion_factor, tolerance).`
///
/// **Note:** The first pmsim's file with index 0 is ignored.
///
/// **Warning:** This function only works with Solid problems with Ux, Uy, and Uz DOFs.
pub fn compare_results(
    mesh: &Mesh,
    schema: &Schema,
    config: &Config,
    dir: &str,
    fn_stem: &str,
    ref_type: ReferenceDataType,
    ref_path: &str,
    tol_displacement: f64,
    tol_stress: f64,
    verbose: usize,
    eps_bar_p: Option<(usize, f64, f64)>,
) -> Result<bool, StrError> {
    // constants
    let dofs = [Dof::Ux, Dof::Uy, Dof::Uz];
    let ndim = mesh.ndim;
    let tensor_vec_dim = 2 * ndim;
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();
    if npoint < 1 {
        return Err("there must be at least one point in the mesh");
    }
    if ncell < 1 {
        return Err("there must be at least one cell");
    }

    // load reference results
    let dat = ReferenceData::load(ref_type, ref_path)?;
    if npoint != dat.actual.npoint() {
        return Err("the number of points in the mesh must equal the corresponding number in the reference data");
    }
    if ncell != dat.actual.ncell() {
        return Err("the number of elements in the mesh must be equal to the reference number of elements (stresses)");
    }

    // stats
    let mut diff_displacement_max = f64::MIN;
    let mut diff_stress_max = f64::MIN;
    let mut diff_eps_bar_p_max = f64::MIN;

    // compare results
    let mut all_good = true;
    let mut elastic_flags_ok = true;
    let (pp, _) = PostProc::new(dir, fn_stem)?;
    if pp.nfile() != dat.actual.nstep() + 1 {
        return Err("the number of steps must equal the reference's number of steps + 1");
    }
    for index in 1..pp.nfile() {
        // set the number of steps in the reference data (where the initial state is absent)
        let step = index - 1;

        // load state
        let fem_state = FemState::read_json(&format!("{}/{}-{}.json", config.out_dir, config.out_fn_stem, index))?;

        if verbose > 0 {
            println!(
                "\nSTEP # {} ===============================================================",
                index
            );
        }

        // check displacements
        if verbose > 0 {
            println!("DISPLACEMENTS");
        }
        for p in 0..npoint {
            for i in 0..ndim {
                let d = schema.dof_number(p, dofs[i])?;
                let a = fem_state.uu[d];
                let b = dat.actual.displacement(step, p, i);
                let (fail, diff) = query_failed(a, b, tol_displacement, verbose);
                diff_displacement_max = f64::max(diff_displacement_max, diff);
                if fail {
                    all_good = false;
                }
            }
            if verbose > 0 {
                println!();
            }
        }

        // check stresses
        if verbose > 0 {
            println!("STRESSES");
        }
        for e in 0..ncell {
            let ngauss = dat.actual.ngauss(step, e);
            if ngauss < 1 {
                return Err("there must be at least on integration point in reference data (stress)");
            }
            let secondary_values = &fem_state.gauss[e];
            for ip in 0..ngauss {
                let local_state = &secondary_values.solid[ip];
                for i in 0..tensor_vec_dim {
                    let a = local_state.stress.vector()[i];
                    let b = if i > 2 {
                        dat.actual.stresses(step, e, ip, i) * SQRT_2 // convert to Mandel
                    } else {
                        dat.actual.stresses(step, e, ip, i)
                    };
                    let (fail, diff) = query_failed(a, b, tol_stress, verbose);
                    diff_stress_max = f64::max(diff_stress_max, diff);
                    if fail {
                        all_good = false;
                    }
                }
                if verbose > 0 {
                    println!();
                }
            }
        }

        // check elastic flags
        if verbose > 0 {
            println!("ELASTIC FLAGS");
        }
        let mut n_elastic = 0;
        for e in 0..ncell {
            let ngauss = dat.actual.ngauss(step, e);
            if ngauss < 1 {
                return Err("there must be at least on integration point in reference data (plast_apex_epbar)");
            }
            let secondary_values = &fem_state.gauss[e];
            for ip in 0..ngauss {
                let local_state = &secondary_values.solid[ip];
                let elastic = dat.actual.elastic(step, e, ip);
                let fail = query_failed_bool(local_state.elastic, elastic, verbose);
                if fail {
                    all_good = false;
                    elastic_flags_ok = false;
                }
                if elastic {
                    n_elastic += 1;
                }
            }
            if verbose > 0 {
                println!();
            }
        }
        if verbose > 0 {
            println!("num elastic = {}", n_elastic);
        }

        // check accumulated plastic strain (eps_bar_p) if requested (von Mise model only)
        if let Some((index, conversion_factor, tolerance)) = eps_bar_p {
            if verbose > 0 {
                println!("ACCUMULATED PLASTIC STRAIN (eps_bar_p)");
            }
            for e in 0..ncell {
                let ngauss = dat.actual.ngauss(step, e);
                if ngauss < 1 {
                    return Err("there must be at least on integration point in reference data (plast_apex_epbar)");
                }
                let secondary_values = &fem_state.gauss[e];
                for ip in 0..ngauss {
                    let local_state = &secondary_values.solid[ip];
                    if local_state.z_set.dim() != NZ_VON_MISES {
                        return Err("the number of internal variables in the local state must equal NZ_VON_MISES");
                    }
                    let a = local_state.z_set[index];
                    let b = dat.actual.eps_bar_p(step, e, ip) * conversion_factor;
                    let (fail, diff) = query_failed(a, b, tolerance, verbose);
                    diff_eps_bar_p_max = f64::max(diff_eps_bar_p_max, diff);
                    if fail {
                        all_good = false;
                    }
                    if verbose > 0 {
                        println!();
                    }
                }
            }
        }
    }
    let s_ok = if elastic_flags_ok { "yes" } else { "no" };
    println!("\ndiff_displacement_max = {:9.2e}", diff_displacement_max);
    println!("diff_stress_max       = {:9.2e}", diff_stress_max);
    if eps_bar_p.is_some() {
        println!("diff_eps_bar_p_max    = {:9.2e}", diff_eps_bar_p_max);
    }
    println!("are elastic flags ok  ? {:>9}\n", s_ok);
    Ok(all_good)
}
