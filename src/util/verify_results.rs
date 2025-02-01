use crate::base::Dof;
use crate::fem::{FemBase, FemState, FileIo};
use crate::util::ReferenceDataSet;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_tensor::SQRT_2;

/// Calculates the error and returns true (fail) if the difference is greater than the tolerance
fn failed(a: f64, b: f64, tol: f64, verbose: bool) -> bool {
    let diff = f64::abs(a - b);
    let fail = diff > tol;
    if verbose {
        let mrk = if fail { "❌" } else { "➖" };
        print!("{:15.6e}{} ", diff, mrk);
    }
    fail
}

/// Compares the FEM results (displacement, stress, strain) against reference data
///
/// # Input
///
/// * `mesh` -- The mesh
/// * `file_io` -- The file output struct
/// * `ref_filename` -- The filename of the file with the reference results (located in `data/results`)
/// * `tol_displacement` -- A tolerance to compare displacements
/// * `tol_stress` -- A tolerance to compare stresses
/// * `verbose` -- Enables the verbose mode
///
/// **Warning:** This function only works with Solid problems with Ux, Uy, and Uz DOFs.
pub fn verify_results(
    mesh: &Mesh,
    base: &FemBase,
    file_io: &FileIo,
    ref_filename: &str,
    tol_displacement: f64,
    tol_stress: f64,
    verbose: bool,
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
    let reference = ReferenceDataSet::read_json(format!("data/results/{}", ref_filename).as_str())?;

    // compare results
    let mut all_good = true;
    let summary = FileIo::read_json(&file_io.path_summary())?;
    for step in &summary.indices {
        if *step >= reference.all.len() {
            return Err("the number of load steps must match the reference data");
        }
        let compare = &reference.all[*step];
        if npoint != compare.displacement.len() {
            return Err("the number of points in the mesh must equal the corresponding number in the reference data");
        }
        if ncell != compare.stresses.len() {
            return Err(
                "the number of elements in the mesh must be equal to the reference number of elements (stresses)",
            );
        }

        // load state
        let fem_state = FemState::read_json(&file_io.path_state(*step))?;

        if verbose {
            println!(
                "\nSTEP # {} ===============================================================",
                step
            );
        }

        // check displacements
        if verbose {
            println!("ERROR ON DISPLACEMENTS");
        }
        for p in 0..npoint {
            if ndim != compare.displacement[p].len() {
                return Err("the space dimension (ndim) must equal the reference data ndim");
            }
            for i in 0..ndim {
                let eq = base.equations.eq(p, dofs[i]).unwrap();
                let a = fem_state.uu[eq];
                let b = compare.displacement[p][i];
                if failed(a, b, tol_displacement, verbose) {
                    all_good = false;
                }
            }
            if verbose {
                println!();
            }
        }

        // check stresses
        println!("ERROR ON STRESSES");
        for e in 0..ncell {
            let ngauss = compare.stresses[e].len();
            if ngauss < 1 {
                return Err("there must be at least on integration point in reference data (stress)");
            }
            let ref_tensor_vec_dim = compare.stresses[e][0].len();
            if tensor_vec_dim != ref_tensor_vec_dim {
                return Err("the number of stress components must equal the reference number of stress components");
            }
            let secondary_values = &fem_state.gauss[e];
            for ip in 0..ngauss {
                let local_state = &secondary_values.solid[ip];
                for i in 0..tensor_vec_dim {
                    let a = local_state.stress.vector()[i];
                    let b = if i > 3 {
                        compare.stresses[e][ip][i] * SQRT_2 // convert to Mandel
                    } else {
                        compare.stresses[e][ip][i]
                    };
                    if failed(a, b, tol_stress, verbose) {
                        all_good = false;
                    }
                }
                println!();
            }
        }
    }
    Ok(all_good)
}
