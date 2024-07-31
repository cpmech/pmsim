use crate::base::DEFAULT_OUT_DIR;
use crate::fem::{FemOutput, FemOutputSummary, FemState};
use crate::material::StressStrainState;
use crate::util::ReferenceDataSet;
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::approx_eq;
use russell_tensor::{Mandel, SQRT_2};

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
/// * `(states, states_ref)` -- The state points at selected `(CellId, IntegrationPointId)`
pub fn check_displacements_and_stresses(
    mesh: &Mesh,
    name: &str,
    ref_filename: &str,
    extract: (usize, usize),
    tol_displacement: f64,
    tol_stress: f64,
) -> Result<(Vec<StressStrainState>, Vec<StressStrainState>), StrError> {
    // constants
    let ndim = mesh.ndim;
    let ncp = 2 * ndim;
    let mandel = Mandel::new(ncp);
    let npoint = mesh.points.len();
    let ncell = mesh.cells.len();
    if npoint < 1 {
        return Err("there must be at least one point");
    }
    if ncell < 1 {
        return Err("there must be at least one cell");
    }

    // load reference results
    let reference = ReferenceDataSet::read_json(format!("data/results/{}", ref_filename).as_str())?;

    // output arrays
    let mut states = Vec::new();
    let mut states_ref = Vec::new();

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
            return Err(
                "the number of elements in the mesh must be equal to the reference number of elements (stresses)",
            );
        }
        if ncell != compare.strains.len() {
            return Err(
                "the number of elements in the mesh must be equal to the reference number of elements (strains)",
            );
        }

        // load state
        let fem_state = FemState::read_json(&FemOutput::path_state(DEFAULT_OUT_DIR, name, *step))?;

        println!("\n===================================================================================");
        println!("STEP # {}", step);

        // check displacements
        println!("ERROR ON DISPLACEMENTS");
        for p in 0..npoint {
            if ndim != compare.displacement[p].len() {
                return Err("the space dimension (ndim) must equal the reference data ndim");
            }
            for i in 0..ndim {
                let a = fem_state.uu[ndim * p + i];
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
                return Err("there must be at least on integration point in reference data (stress)");
            }
            let ncp_ref = compare.stresses[e][0].len();
            if ncp != ncp_ref {
                return Err("the number of stress components must equal the reference number of stress components");
            }
            for ip in 0..n_integ_point {
                // extract state at integration point
                let state = fem_state.extract_stresses_and_strains(e, ip).unwrap();

                // check stresses
                for i in 0..ncp {
                    let a = state.stress.vector()[i];
                    let b = if i > 3 {
                        compare.stresses[e][ip][i] * SQRT_2 // convert to Mandel
                    } else {
                        compare.stresses[e][ip][i]
                    };
                    print!("{:13.6e} ", f64::abs(a - b));
                    approx_eq(a, b, tol_stress);
                }
                println!();
            }
        }

        // check strains
        println!("ERROR ON STRAINS");
        for e in 0..ncell {
            let n_integ_point = compare.stresses[e].len();
            let n_integ_point_eps = compare.strains[e].len();
            if n_integ_point_eps != n_integ_point {
                return Err("strain data must have the same number of integration points as stress data");
            }
            let ncp_ref_eps = compare.strains[e][0].len();
            if ncp_ref_eps != ncp {
                return Err("strain data must have the same number of components as stress data");
            }
            for ip in 0..n_integ_point {
                // extract state at integration point
                let state = fem_state.extract_stresses_and_strains(e, ip).unwrap();

                // check strains
                for i in 0..ncp {
                    let a = state.strain().vector()[i];
                    let b = if i > 3 {
                        compare.strains[e][ip][i] * SQRT_2 // convert to Mandel
                    } else {
                        compare.strains[e][ip][i]
                    };
                    // TODO: check this
                    print!("{:13.6e}({:13.6e}) ", a, b);
                    // print!("{:13.6e} ", f64::abs(a - b));
                    // approx_eq(a, b, tol_stress);
                }
                println!();
            }
        }

        // save the results at selected integration points
        for e in 0..ncell {
            let n_integ_point = compare.stresses[e].len();
            let with_optional = true;
            for ip in 0..n_integ_point {
                if e == extract.0 && ip == extract.1 {
                    let state = fem_state.extract_stresses_and_strains(e, ip).unwrap();
                    let mut state_ref = StressStrainState::new(mandel, 0, with_optional);
                    for i in 0..ncp {
                        if i > 3 {
                            state_ref.stress.vector_mut()[i] = compare.stresses[e][ip][i] * SQRT_2;
                            state_ref.strain_mut().vector_mut()[i] = compare.strains[e][ip][i] * SQRT_2;
                        } else {
                            state_ref.stress.vector_mut()[i] = compare.stresses[e][ip][i];
                            state_ref.strain_mut().vector_mut()[i] = compare.strains[e][ip][i];
                        }
                    }
                    states.push(state.clone());
                    states_ref.push(state_ref);
                }
            }
        }
    }
    Ok((states, states_ref))
}
