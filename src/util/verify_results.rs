use crate::base::DEFAULT_OUT_DIR;
use crate::fem::{FemOutput, FemOutputSummary, FemState};
use crate::util::ReferenceDataSet;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::approx_eq;

/// Validates the FEM results (displacement, stress, strain) against reference data
///
/// # Input
///
/// * `mesh` -- The mesh
/// * `name` -- The FemOutput file name
/// * `ref_filename` -- The reference results filename (in `data/results`)
/// * `tol_displacement` -- A tolerance to compare displacements
/// * `tol_stress` -- A tolerance to compare stresses
///
/// # Output
///
/// * `(states, states_ref)` -- The state points at selected `(CellId, IntegrationPointId)`
pub fn validate_results(
    mesh: &Mesh,
    name: &str,
    ref_filename: &str,
    tol_displacement: f64,
    verbose: bool,
) -> Result<(), StrError> {
    // constants
    let ndim = mesh.ndim;
    let npoint = mesh.points.len();
    if npoint < 1 {
        return Err("there must be at least one point in the mesh");
    }

    // load reference results
    let reference = ReferenceDataSet::read_json(format!("data/results/{}", ref_filename).as_str())?;

    // compare results
    let summary = FemOutputSummary::read_json(&FemOutput::path_summary(DEFAULT_OUT_DIR, name))?;
    for step in &summary.indices {
        if *step >= reference.all.len() {
            return Err("the number of load steps must match the reference data");
        }
        let compare = &reference.all[*step];
        if npoint != compare.displacement.len() {
            return Err("the number of points in the mesh must equal the corresponding number in the reference data");
        }

        // load state
        let fem_state = FemState::read_json(&FemOutput::path_state(DEFAULT_OUT_DIR, name, *step))?;

        if verbose {
            println!("\nSTEP # {} ===================================================", step);
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
                let a = fem_state.uu[ndim * p + i];
                let b = compare.displacement[p][i];
                if verbose {
                    print!("{:13.6e} ", f64::abs(a - b));
                }
                approx_eq(a, b, tol_displacement);
            }
            if verbose {
                println!();
            }
        }
    }
    Ok(())
}
