use crate::base::DEFAULT_OUT_DIR;
use crate::fem::{FemOutput, FemOutputSummary, FemState};
use crate::util::ReferenceDataSet;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::approx_eq;
use russell_tensor::{Tensor2, SQRT_2};

/// Checks displacements and results against reference data
///
/// # Input
///
/// * `mesh` -- The mesh
/// * `name` -- The FemOutput file name
/// * `ref_filename` -- The reference results filename (in `data/results`)
/// * `extract` -- A pair of (CellId, IntegrationPointId)
/// * `tol_displacement` -- A tolerance to compare displacements
/// * `tol_stress` -- A tolerance to compare stresses
///
/// # Output
///
/// * `(stresses, ref_stresses)` -- The stress paths at selected (CellId, IntegrationPointId)
pub fn check_displacements_and_stresses(
    mesh: &Mesh,
    name: &str,
    ref_filename: &str,
    extract: (usize, usize),
    tol_displacement: f64,
    tol_stress: f64,
) -> Result<(Vec<Tensor2>, Vec<Tensor2>), StrError> {
    // constants
    let ndim = mesh.ndim;
    let nsigma = 2 * ndim;
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();
    if npoint < 1 {
        return Err("there must be at least one point");
    }
    if ncell < 1 {
        return Err("there must be at least one cell");
    }

    // output arrays
    let mut stresses = Vec::new();
    let mut ref_stresses = Vec::new();

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
        if ncell != compare.stresses.len() {
            return Err("the number of elements in the mesh must be equal to the reference number of elements");
        }

        // load state
        let state = FemState::read_json(&FemOutput::path_state(DEFAULT_OUT_DIR, name, *step))?;

        println!("\n===================================================================================");
        println!("STEP # {}", step);

        // check displacements
        println!("ERROR ON DISPLACEMENTS");
        for p in 0..npoint {
            if ndim != compare.displacement[p].len() {
                return Err("the space dimension (ndim) must equal the reference data ndim");
            }
            for i in 0..ndim {
                let a = state.uu[ndim * p + i];
                let b = compare.displacement[p][i];
                print!("{:13.6e} ", f64::abs(a - b));
                approx_eq(a, b, tol_displacement);
            }
            println!();
        }

        // check stresses
        println!("ERROR ON STRESSES");
        for e in 0..ncell {
            let n_integ_point = compare.stresses[e].len();
            if n_integ_point < 1 {
                return Err("there must be at least on integration point in reference data");
            }
            let ref_nsigma = compare.stresses[e][0].len();
            if nsigma != ref_nsigma {
                return Err("the number of stress components must equal the reference number of stress components");
            }
            for ip in 0..n_integ_point {
                // check sigma
                let state = state.extract_stresses_and_strains(e, ip).unwrap();
                for i in 0..nsigma {
                    let a = state.sigma.vector()[i];
                    let b = if i > 3 {
                        compare.stresses[e][ip][i] * SQRT_2
                    } else {
                        compare.stresses[e][ip][i]
                    };
                    print!("{:13.6e} ", f64::abs(a - b));
                    approx_eq(a, b, tol_stress);
                }
                println!();

                // save stress path at selected (CellId, IntegrationPointId)
                if e == extract.0 && ip == extract.1 {
                    stresses.push(state.sigma.clone());
                    let mut sigma = Tensor2::new_sym(true);
                    for i in 0..nsigma {
                        if i > 3 {
                            sigma.vector_mut()[i] = compare.stresses[e][ip][i] * SQRT_2;
                        } else {
                            sigma.vector_mut()[i] = compare.stresses[e][ip][i];
                        }
                    }
                    ref_stresses.push(sigma);
                }
            }
        }
    }
    Ok((stresses, ref_stresses))
}
