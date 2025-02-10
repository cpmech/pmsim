use crate::base::Dof;
use crate::fem::{FemBase, FemState, FileIo};
use crate::util::ReferenceDataSet;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_tensor::SQRT_2;

pub(crate) trait ReferenceDataTrait {
    /// Returns the displacements
    ///
    /// Size: `[npoint][ndim]`
    fn displacement(&self) -> &Vec<Vec<f64>>;

    /// Returns the stresses (standard components)
    ///
    /// Size: `[nele][ngauss][n_components]`
    fn stresses(&self) -> &Vec<Vec<Vec<f64>>>;
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

/// Compares the FEM results (displacement, stress, strain) against reference data
///
/// # Input
///
/// * `mesh` -- The mesh
/// * `file_io` -- The file output struct
/// * `ref_full_path` -- The full path of the file with the reference results
/// * `tol_displacement` -- A tolerance to compare displacements
/// * `tol_stress` -- A tolerance to compare stresses
/// * `verbose` -- Enables the verbose mode:
///   - 0 => no output
///   - 1 => shows error
///   - 2 => shows values and error
///
/// **Note:** The first pmsim's file with index 0 is ignored.
///
/// **Warning:** This function only works with Solid problems with Ux, Uy, and Uz DOFs.
pub fn compare_results(
    mesh: &Mesh,
    base: &FemBase,
    file_io: &FileIo,
    ref_full_path: &str,
    tol_displacement: f64,
    tol_stress: f64,
    verbose: usize,
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
    let reference = ReferenceDataSet::read_json(ref_full_path)?;

    // stats
    let mut diff_displacement_max = f64::MIN;
    let mut diff_stress_max = f64::MIN;

    // compare results
    let mut all_good = true;
    let summary = FileIo::read_json(&file_io.path_summary())?;
    if summary.indices.len() != reference.all.len() + 1 {
        return Err("the number of steps must equal the reference's number of steps + 1");
    }
    for step in 1..summary.indices.len() {
        let compare = &reference.all[step - 1];
        if npoint != compare.displacement().len() {
            return Err("the number of points in the mesh must equal the corresponding number in the reference data");
        }
        if ncell != compare.stresses().len() {
            return Err(
                "the number of elements in the mesh must be equal to the reference number of elements (stresses)",
            );
        }

        // load state
        let fem_state = FemState::read_json(&file_io.path_state(step))?;

        if verbose > 0 {
            println!(
                "\nSTEP # {} ===============================================================",
                step
            );
        }

        // check displacements
        if verbose > 0 {
            println!("ERROR ON DISPLACEMENTS");
        }
        for p in 0..npoint {
            if ndim != compare.displacement()[p].len() {
                return Err("the space dimension (ndim) must equal the reference data ndim");
            }
            for i in 0..ndim {
                let eq = base.equations.eq(p, dofs[i]).unwrap();
                let a = fem_state.uu[eq];
                let b = compare.displacement()[p][i];
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
            println!("ERROR ON STRESSES");
        }
        for e in 0..ncell {
            let ngauss = compare.stresses()[e].len();
            if ngauss < 1 {
                return Err("there must be at least on integration point in reference data (stress)");
            }
            let ref_tensor_vec_dim = compare.stresses()[e][0].len();
            if tensor_vec_dim != ref_tensor_vec_dim {
                return Err("the number of stress components must equal the reference number of stress components");
            }
            let secondary_values = &fem_state.gauss[e];
            for ip in 0..ngauss {
                let local_state = &secondary_values.solid[ip];
                for i in 0..tensor_vec_dim {
                    let a = local_state.stress.vector()[i];
                    let b = if i > 2 {
                        compare.stresses()[e][ip][i] * SQRT_2 // convert to Mandel
                    } else {
                        compare.stresses()[e][ip][i]
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
    }
    if verbose > 0 {
        println!("\ndiff_displacement_max = {:9.2e}", diff_displacement_max);
        println!("diff_stress_max       = {:9.2e}", diff_stress_max);
        println!();
    }
    Ok(all_good)
}
