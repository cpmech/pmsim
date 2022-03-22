use super::{Dof, EquationNumbers, State};
use crate::elements::Element;
use crate::StrError;
use russell_lab::{mat_max_abs_diff, Matrix};
use russell_tensor::SQRT_2;
use serde::Deserialize;
use serde_json;
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

const VALIDATOR_DEFAULT_TOL_K_MATRIX: f64 = 1e-12;
const VALIDATOR_DEFAULT_TOL_DISPLACEMENT: f64 = 1e-12;
const VALIDATOR_DEFAULT_TOL_STRESS: f64 = 1e-12;

/// Holds results from iterations
#[derive(Clone, Debug, Deserialize)]
pub struct ValidatorIteration {
    /// Iteration number
    #[serde(rename(deserialize = "it"))]
    pub iteration: usize,

    /// Relative residual
    #[serde(rename(deserialize = "resrel"))]
    pub relative_residual: f64,

    /// Absolute residual
    #[serde(rename(deserialize = "resid"))]
    pub absolute_residual: f64,
}

/// Holds numerical results
#[derive(Clone, Debug, Deserialize)]
pub struct ValidatorResults {
    /// All stiffness matrices (nele,nu,nu)
    #[serde(default)]
    #[serde(rename(deserialize = "Kmats"))]
    pub kk_matrices: Vec<Matrix>,

    /// Displacements at nodes (npoint,ndim)
    #[serde(default)]
    #[serde(rename(deserialize = "disp"))]
    pub displacements: Vec<Vec<f64>>,

    /// All stresses @ all ips (nele,nip,nsigma)
    ///
    /// NOTE: these are "real" components and the order is (sx, sy, sxy, {sz}), different than Tensor2
    #[serde(default)]
    pub stresses: Vec<Vec<Vec<f64>>>,

    /// Load factor
    #[serde(default)]
    #[serde(rename(deserialize = "loadfactor"))]
    pub load_factor: f64,

    /// Iterations data
    #[serde(default)]
    pub iterations: Vec<ValidatorIteration>,
}

/// Holds results for comparisons (checking/tests)
#[derive(Clone, Debug, Deserialize)]
pub struct Validator {
    /// Results from output-, time-, or load-steps
    pub steps: Vec<ValidatorResults>,

    /// Tolerance to compare K matrix components
    #[serde(skip)]
    pub tol_kk_matrix: f64,

    /// Tolerance to compare displacement components
    #[serde(skip)]
    pub tol_displacement: f64,

    /// Tolerance to compare stress components
    #[serde(skip)]
    pub tol_stress: f64,
}

impl Validator {
    /// Returns a new Validator from JSON string
    pub fn from_str(cmp: &str) -> Result<Self, StrError> {
        let mut res: Validator = serde_json::from_str(&cmp).map_err(|op| {
            println!("ERROR: {}", op);
            return "serde_json failed";
        })?;
        res.tol_kk_matrix = VALIDATOR_DEFAULT_TOL_K_MATRIX;
        res.tol_displacement = VALIDATOR_DEFAULT_TOL_DISPLACEMENT;
        res.tol_stress = VALIDATOR_DEFAULT_TOL_STRESS;
        Ok(res)
    }

    /// Reads a JSON file containing the comparison results
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let file = File::open(&path).map_err(|_| "file not found")?;
        let reader = BufReader::new(file);
        let mut res: Validator = serde_json::from_reader(reader).map_err(|op| {
            println!("ERROR: {}", op);
            return "serde_json failed";
        })?;
        res.tol_kk_matrix = VALIDATOR_DEFAULT_TOL_K_MATRIX;
        res.tol_displacement = VALIDATOR_DEFAULT_TOL_DISPLACEMENT;
        res.tol_stress = VALIDATOR_DEFAULT_TOL_STRESS;
        Ok(res)
    }

    /// Compare K matrices against results at a fixed time- or load-step
    ///
    /// Returns "OK" if all values are approximately equal under tol_kk_matrix.
    pub fn compare_kk_matrices(&self, step: usize, state: &State, elements: &mut Vec<Element>) -> String {
        if step >= self.steps.len() {
            return "reference results for the step are not available".to_string();
        }
        let cmp = &self.steps[step];
        let nele = elements.len();
        for element_id in 0..nele {
            if element_id >= cmp.kk_matrices.len() {
                return format!("element {}: reference K matrix is not available", element_id);
            }
            let element = &mut elements[element_id];
            if let Err(e) = element.base.calc_local_kk_matrix(&state.elements[element_id], true) {
                return format!("element {}: calc_local_kk_matrix failed: {}", element_id, e);
            }
            let kk = element.base.get_local_kk_matrix();
            let reference = &cmp.kk_matrices[element_id];
            match mat_max_abs_diff(&kk, &reference) {
                Ok((i, j, max_abs_diff)) => {
                    if max_abs_diff > self.tol_kk_matrix {
                        return format!(
                            "element {}: K{}{} component is greater than tolerance. max_abs_diff = {:e}",
                            element_id, i, j, max_abs_diff
                        );
                    }
                }
                Err(e) => return format!("element {}: mat_max_abs_diff failed: {}", element_id, e),
            }
        }
        "OK".to_string()
    }

    /// Compares displacements against results at a fixed time- or load-step
    ///
    /// Returns "OK" if all values are approximately equal under tol_displacement.
    pub fn compare_displacements(
        &self,
        step: usize,
        state: &State,
        equations: &EquationNumbers,
        two_dim: bool,
    ) -> String {
        if step >= self.steps.len() {
            return "reference results for the step are not available".to_string();
        }
        let cmp = &self.steps[step];
        let n_disp_components = if two_dim { 2 } else { 3 };
        let npoint = equations.npoint();
        for point_id in 0..npoint {
            if point_id >= cmp.displacements.len() {
                return format!("point {}: reference displacement is not available", point_id);
            }
            let reference = &cmp.displacements[point_id];
            if reference.len() != n_disp_components {
                return format!(
                    "point {}: reference displacement has incompatible number of components. {}(wrong) != {}",
                    point_id,
                    reference.len(),
                    n_disp_components
                );
            }
            let ux = match equations.number(point_id, Dof::Ux) {
                Some(eq) => {
                    if eq >= state.system_xx.dim() {
                        return format!(
                            "point {}: state does not have equation {} corresponding to ux",
                            point_id, eq
                        );
                    }
                    state.system_xx[eq]
                }
                None => return format!("point {}: state does not have ux", point_id),
            };
            let uy = match equations.number(point_id, Dof::Uy) {
                Some(eq) => {
                    if eq >= state.system_xx.dim() {
                        return format!(
                            "point {}: state does not have equation {} corresponding to uy",
                            point_id, eq
                        );
                    }
                    state.system_xx[eq]
                }
                None => return format!("point {}: state does not have uy", point_id),
            };
            let diff_ux = f64::abs(ux - reference[0]);
            if diff_ux > self.tol_displacement {
                return format!(
                    "point {}: ux is greater than tolerance. |ux - reference| = {:e}",
                    point_id, diff_ux
                );
            }
            let diff_uy = f64::abs(uy - reference[1]);
            if diff_uy > self.tol_displacement {
                return format!(
                    "point {}: uy is greater than tolerance. |uy - reference| = {:e}",
                    point_id, diff_uy
                );
            }
            if !two_dim {
                let uz = match equations.number(point_id, Dof::Uz) {
                    Some(eq) => {
                        if eq >= state.system_xx.dim() {
                            return format!(
                                "point {}: state does not have equation {} corresponding to uz",
                                point_id, eq
                            );
                        }
                        state.system_xx[eq]
                    }
                    None => return format!("point {}: state does not have uz", point_id),
                };
                let diff_uz = f64::abs(uz - reference[2]);
                if diff_uz > self.tol_displacement {
                    return format!(
                        "point {}: uz is greater than tolerance. |uz - reference| = {:e}",
                        point_id, diff_uz
                    );
                }
            }
        }
        "OK".to_string()
    }

    /// Compares stresses against results at a fixed time- or load-step
    ///
    /// Returns "OK" if all values are approximately equal under tol_stress.
    pub fn compare_stresses(&self, step: usize, state: &State, two_dim: bool) -> String {
        if step >= self.steps.len() {
            return "reference results for the step are not available".to_string();
        }
        let cmp = &self.steps[step];
        let n_stress_components = if two_dim { 3 } else { 4 };
        let nele = state.elements.len();
        for element_id in 0..nele {
            let element = &state.elements[element_id];
            let n_integ_point = element.stress.len();
            if n_integ_point < 1 {
                return format!("element {}: state does not have integration point stress", element_id);
            }
            if element_id >= cmp.stresses.len() {
                return format!("element {}: reference stress is not available", element_id);
            }
            for index_ip in 0..n_integ_point {
                if index_ip >= cmp.stresses[element_id].len() {
                    return format!(
                        "element {}: integration point {}: reference stress is not available",
                        element_id, index_ip
                    );
                }
                let reference = &cmp.stresses[element_id][index_ip];
                if reference.len() != n_stress_components {
                    return format!(
                            "element {}: integration point {}: reference stress has incompatible number of components. {}(wrong) != {}",
                            element_id, index_ip, reference.len(), n_stress_components,
                        );
                }
                let sigma = &element.stress[index_ip].sigma;
                let sx = sigma.vec[0];
                let sy = sigma.vec[1];
                let sxy = sigma.vec[3] / SQRT_2;
                let diff_sx = f64::abs(sx - reference[0]);
                if diff_sx > self.tol_stress {
                    return format!(
                        "element {}: integration point {}: sx is greater than tolerance. |sx - reference| = {:e}",
                        element_id, index_ip, diff_sx
                    );
                }
                let diff_sy = f64::abs(sy - reference[1]);
                if diff_sy > self.tol_stress {
                    return format!(
                        "element {}: integration point {}: sy is greater than tolerance. |sy - reference| = {:e}",
                        element_id, index_ip, diff_sy
                    );
                }
                let diff_sxy = f64::abs(sxy - reference[2]);
                if diff_sxy > self.tol_stress {
                    return format!(
                        "element {}: integration point {}: sxy is greater than tolerance. |sxy - reference| = {:e}",
                        element_id, index_ip, diff_sxy
                    );
                }
                if !two_dim {
                    let sz = sigma.vec[2];
                    let diff_sz = f64::abs(sz - reference[3]);
                    if diff_sz > self.tol_stress {
                        return format!(
                            "element {}: integration point {}: sz is greater than tolerance. |sz - reference| = {:e}",
                            element_id, index_ip, diff_sz
                        );
                    }
                }
            }
        }
        "OK".to_string()
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Validator, ValidatorIteration};
    use crate::elements::Element;
    use crate::simulation::{
        Configuration, Dof, ElementConfig, EquationNumbers, Initializer, ParamSolid, ParamStressStrain, State,
        StateElement, StateStress,
    };
    use crate::StrError;
    use gemlab::mesh::Mesh;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;
    use russell_tensor::Tensor2;

    #[test]
    fn serialize_handles_errors() {
        assert_eq!(Validator::from_str("").err(), Some("serde_json failed"));
        assert_eq!(Validator::read_json("").err(), Some("file not found"));
        assert_eq!(
            Validator::read_json("./data/validation/z_wrong_data_for_test.json").err(),
            Some("serde_json failed")
        );
    }

    #[test]
    fn clone_debug_and_serialize_work() -> Result<(), StrError> {
        let ite = ValidatorIteration {
            iteration: 1,
            relative_residual: 2.0,
            absolute_residual: 3.0,
        };
        let cloned = ite.clone();
        assert_eq!(cloned.iteration, 1);
        let ite: ValidatorIteration = serde_json::from_str(r#"{"it":1,"resrel":2,"resid":3}"#).unwrap();
        assert_eq!(ite.iteration, 1);
        assert_eq!(
            format!("{:?}", ite),
            "ValidatorIteration { iteration: 1, relative_residual: 2.0, absolute_residual: 3.0 }"
        );

        let val = Validator::from_str(
            r#"{ "steps":
          [
            {
              "Kmats": [
                { "nrow":2, "ncol":2, "data":[
                  1, 2,
                  3, 4
                ] },
                { "nrow":2, "ncol":2, "data":[
                  10, 20,
                  30, 40
                ] }
              ],
              "disp": [
                [11, 21],
                [12, 22]
              ],
              "stresses": [
                [
                  [100, 101, 102, 103]
                ],
                [
                  [200, 201, 202, 203]
                ]
              ]
            }
          ]
        }"#,
        )?;
        assert_eq!(val.steps.len(), 1);
        let res = &val.steps[0];
        assert_eq!(res.kk_matrices.len(), 2);
        assert_eq!(
            format!("{}", res.kk_matrices[0]),
            "┌     ┐\n\
             │ 1 2 │\n\
             │ 3 4 │\n\
             └     ┘"
        );
        assert_eq!(
            format!("{}", res.kk_matrices[1]),
            "┌       ┐\n\
             │ 10 20 │\n\
             │ 30 40 │\n\
             └       ┘"
        );
        assert_eq!(res.displacements.len(), 2);
        assert_vec_approx_eq!(res.displacements[0], [11.0, 21.0], 1e-15);
        assert_vec_approx_eq!(res.displacements[1], [12.0, 22.0], 1e-15);
        assert_eq!(res.stresses.len(), 2);
        assert_vec_approx_eq!(res.stresses[0][0], [100.0, 101.0, 102.0, 103.0], 1e-15);
        assert_vec_approx_eq!(res.stresses[1][0], [200.0, 201.0, 202.0, 203.0], 1e-15);

        let cloned = val.clone();
        assert_eq!(cloned.steps.len(), 1);
        let c_res = &cloned.steps[0];
        assert_eq!(c_res.kk_matrices.len(), 2);
        assert_eq!(
            format!("{}", c_res.kk_matrices[0]),
            "┌     ┐\n\
             │ 1 2 │\n\
             │ 3 4 │\n\
             └     ┘"
        );
        assert_eq!(
            format!("{}", c_res.kk_matrices[1]),
            "┌       ┐\n\
             │ 10 20 │\n\
             │ 30 40 │\n\
             └       ┘"
        );
        assert_eq!(c_res.displacements.len(), 2);
        assert_vec_approx_eq!(c_res.displacements[0], [11.0, 21.0], 1e-15);
        assert_vec_approx_eq!(c_res.displacements[1], [12.0, 22.0], 1e-15);
        assert_eq!(c_res.stresses.len(), 2);
        assert_vec_approx_eq!(c_res.stresses[0][0], [100.0, 101.0, 102.0, 103.0], 1e-15);
        assert_vec_approx_eq!(c_res.stresses[1][0], [200.0, 201.0, 202.0, 203.0], 1e-15);

        assert_eq!(format!("{:?}", val), "Validator { steps: [ValidatorResults { kk_matrices: [NumMatrix { nrow: 2, ncol: 2, data: [1.0, 2.0, 3.0, 4.0] }, NumMatrix { nrow: 2, ncol: 2, data: [10.0, 20.0, 30.0, 40.0] }], displacements: [[11.0, 21.0], [12.0, 22.0]], stresses: [[[100.0, 101.0, 102.0, 103.0]], [[200.0, 201.0, 202.0, 203.0]]], load_factor: 0.0, iterations: [] }], tol_kk_matrix: 1e-12, tol_displacement: 1e-12, tol_stress: 1e-12 }");
        Ok(())
    }

    #[test]
    fn from_json_works() -> Result<(), StrError> {
        let val = Validator::from_str(
            r#"{ "steps":
          [
            {
              "loadfactor" : 1,
              "iterations" : [
                { "it":1, "resrel": 1, "resid": 2}
              ],
              "disp" : [
                [1, 2],
                [1, 2]
              ],
              "stresses" : [
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ],
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ]
              ]
            },
            {
              "loadfactor" : 2,
              "iterations" : [
                { "it":1, "resrel": 1, "resid": 2},
                { "it":2, "resrel": 1, "resid": 2}
              ],
              "disp" : [
                [1, 2],
                [1, 2]
              ],
              "stresses" : [
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ],
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ]
              ]
            },
            {
              "loadfactor" : 3,
              "iterations" : [
                { "it":1, "resrel": 1, "resid": 2},
                { "it":2, "resrel": 1, "resid": 2},
                { "it":3, "resrel": 1, "resid": 2}
              ],
              "disp" : [
                [1, 2],
                [1, 2]
              ],
              "stresses" : [
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ],
                [
                  [1, 2, 3, 4],
                  [1, 2, 3, 4]
                ]
              ]
            }
          ]
        }"#,
        )?;
        assert_eq!(val.steps.len(), 3);
        assert_eq!(val.steps[0].iterations.len(), 1);
        assert_eq!(val.steps[1].iterations.len(), 2);
        assert_eq!(val.steps[2].iterations.len(), 3);
        Ok(())
    }

    #[test]
    fn read_json_works() -> Result<(), StrError> {
        let val = Validator::read_json("./data/validation/bhatti_1_6.json")?;
        assert_eq!(val.steps.len(), 1);
        assert_eq!(val.steps[0].kk_matrices.len(), 4);
        assert_eq!(val.steps[0].displacements.len(), 6);
        assert_eq!(val.steps[0].stresses.len(), 4);
        Ok(())
    }

    #[test]
    fn compare_kk_matrices_captures_errors() -> Result<(), StrError> {
        /* Example 1.6 from [@bhatti] page 32

         Solid bracket with thickness = 0.25

                      1    load                connectivity:
         y=2.0  fixed *'-,__                    eid : vertices
                      |     '-,_  3   load        0 :  0, 2, 3
         y=1.5 - - -  |        ,'*-,__            1 :  3, 1, 0
                      |  1   ,'  |    '-,_  5     2 :  2, 4, 5
         y=1.0 - - -  |    ,'    |  3   ,-'*      3 :  5, 3, 2
                      |  ,'  0   |   ,-'   |
                      |,'        |,-'   2  |   constraints:
         y=0.0  fixed *----------*---------*     fixed on x and y
                      0          2         4
                     x=0.0     x=2.0     x=4.0

        # References

        [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
                  and Applications, Wiley, 700p.
        */
        let mesh = Mesh::from_text_file("./data/meshes/bhatti_1_6.msh")?;

        let mut config = Configuration::new(&mesh);
        config.plane_stress(true)?.thickness(0.24)?; // << wrong
        let param = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0,
                poisson: 0.2,
            },
        };
        config.elements(1, ElementConfig::Solid(param, None))?;
        let initializer = Initializer::new(&config)?;

        let element_0 = Element::new(&config, 0)?;
        let state = State {
            elements: vec![element_0.base.new_state(&initializer)?],
            system_xx: Vector::new(0),
            system_yy: Vector::new(0),
        };
        let mut elements = vec![element_0];

        let val = Validator::from_str(r#"{ "steps": [] }"#)?;
        assert_eq!(
            val.compare_kk_matrices(0, &state, &mut elements),
            "reference results for the step are not available"
        );

        let val = Validator::from_str(r#"{ "steps": [ { "Kmats":[] } ] }"#)?;
        assert_eq!(
            val.compare_kk_matrices(0, &state, &mut elements),
            "element 0: reference K matrix is not available"
        );

        let val = Validator::from_str(r#"{ "steps": [ { "Kmats":[ {"nrow":2,"ncol":2,"data":[1,2,3,4]} ] } ] }"#)?;
        assert_eq!(
            val.compare_kk_matrices(0, &state, &mut elements),
            "element 0: mat_max_abs_diff failed: matrices are incompatible"
        );

        let val = Validator::read_json("./data/validation/bhatti_1_6.json")?;
        let res = val.compare_kk_matrices(0, &state, &mut elements);
        assert_eq!(
            &res[..66],
            "element 0: K33 component is greater than tolerance. max_abs_diff ="
        );
        Ok(())
    }

    #[test]
    fn compare_kk_matrices_works() -> Result<(), StrError> {
        /* Example 1.6 from [@bhatti] page 32

         Solid bracket with thickness = 0.25

                      1    load                connectivity:
         y=2.0  fixed *'-,__                    eid : vertices
                      |     '-,_  3   load        0 :  0, 2, 3
         y=1.5 - - -  |        ,'*-,__            1 :  3, 1, 0
                      |  1   ,'  |    '-,_  5     2 :  2, 4, 5
         y=1.0 - - -  |    ,'    |  3   ,-'*      3 :  5, 3, 2
                      |  ,'  0   |   ,-'   |
                      |,'        |,-'   2  |   constraints:
         y=0.0  fixed *----------*---------*     fixed on x and y
                      0          2         4
                     x=0.0     x=2.0     x=4.0

        # References

        [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
                  and Applications, Wiley, 700p.
        */
        let mesh = Mesh::from_text_file("./data/meshes/bhatti_1_6.msh")?;

        let mut config = Configuration::new(&mesh);
        config.plane_stress(true)?.thickness(0.25)?;
        let param = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0,
                poisson: 0.2,
            },
        };
        config.elements(1, ElementConfig::Solid(param, None))?;
        let initializer = Initializer::new(&config)?;

        let element_0 = Element::new(&config, 0)?;
        let state = State {
            elements: vec![element_0.base.new_state(&initializer)?],
            system_xx: Vector::new(0),
            system_yy: Vector::new(0),
        };
        let mut elements = vec![element_0];

        let val = Validator::read_json("./data/validation/bhatti_1_6.json")?;
        let res = val.compare_kk_matrices(0, &state, &mut elements);
        assert_eq!(res, "OK");
        Ok(())
    }

    #[test]
    fn compare_displacements_captures_errors_2d() -> Result<(), StrError> {
        let two_dim = true;
        let mut equations = EquationNumbers::new(1);
        let mut state = State {
            elements: Vec::new(),
            system_xx: Vector::new(0), // << incorrect
            system_yy: Vector::new(0), // << incorrect
        };

        let val = Validator::from_str(r#"{ "steps": [] }"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "reference results for the step are not available"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: reference displacement is not available"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: reference displacement has incompatible number of components. 1(wrong) != 2"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have ux"
        );

        equations.activate_equation(0, Dof::Ux);
        state.system_xx = Vector::new(2);

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have uy"
        );

        equations.activate_equation(0, Dof::Uy);
        state.system_xx = Vector::new(0); // << incorrect

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have equation 0 corresponding to ux"
        );

        state.system_xx = Vector::new(1); // << incorrect

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have equation 1 corresponding to uy"
        );

        let neq = equations.nequation();
        let mut state = State {
            elements: Vec::new(),
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
        };
        state.system_xx[0] = 3.0;
        state.system_xx[1] = 4.0;

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: ux is greater than tolerance. |ux - reference| = 2e0"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[3,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: uy is greater than tolerance. |uy - reference| = 2e0"
        );
        Ok(())
    }

    #[test]
    fn compare_displacements_captures_errors_3d() -> Result<(), StrError> {
        let two_dim = false;
        let mut equations = EquationNumbers::new(1);
        let state = State {
            elements: Vec::new(),
            system_xx: Vector::new(2), // << incorrect
            system_yy: Vector::new(2), // << incorrect
        };

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: reference displacement has incompatible number of components. 2(wrong) != 3"
        );

        equations.activate_equation(0, Dof::Ux);
        equations.activate_equation(0, Dof::Uy);

        let val = Validator::from_str(r#"{"steps":[{"disp":[[0,0,0]],"stresses":[[[1,2,3,4]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have uz"
        );

        equations.activate_equation(0, Dof::Uz);

        let val = Validator::from_str(r#"{"steps":[{"disp":[[0,0,0]],"stresses":[[[1,2,3,4]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: state does not have equation 2 corresponding to uz"
        );

        let neq = equations.nequation();
        let mut state = State {
            elements: Vec::new(),
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
        };
        state.system_xx[0] = 3.0;
        state.system_xx[1] = 4.0;
        state.system_xx[2] = 6.0;

        let val = Validator::from_str(r#"{"steps":[{"disp":[[3,4,4]],"stresses":[[[1,2,3,4]]]}]}"#)?;
        assert_eq!(
            val.compare_displacements(0, &state, &equations, two_dim),
            "point 0: uz is greater than tolerance. |uz - reference| = 2e0"
        );
        Ok(())
    }

    #[test]
    fn compare_stresses_handles_missing_data() -> Result<(), StrError> {
        let two_dim = true;
        let equations = EquationNumbers::new(1);
        let neq = equations.nequation();
        let state = State {
            elements: vec![StateElement {
                seepage: Vec::new(),
                stress: Vec::new(),
            }],
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
        };
        let val = Validator::from_str(r#"{"steps":[{"disp":[[0,0]],"stresses":[[[1,2,3]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: state does not have integration point stress"
        );
        Ok(())
    }

    #[test]
    fn compare_stresses_captures_errors_2d() -> Result<(), StrError> {
        let two_dim = true;
        let mut equations = EquationNumbers::new(1);
        equations.activate_equation(0, Dof::Ux);
        equations.activate_equation(0, Dof::Uy);
        let neq = equations.nequation();
        let sigma = StateStress {
            internal_values: Vec::new(),
            sigma: Tensor2::from_matrix(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], true, two_dim)?,
        };
        let mut state = State {
            elements: vec![StateElement {
                seepage: Vec::new(),
                stress: vec![sigma.clone(), sigma.clone()],
            }],
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
        };
        state.system_xx[0] = 1.0;
        state.system_xx[1] = 2.0;

        let val = Validator::from_str(r#"{ "steps": [] }"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "reference results for the step are not available"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: reference stress is not available"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[0,0,0]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 1: reference stress is not available"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[0,0,0],[1,2]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 1: reference stress has incompatible number of components. 2(wrong) != 3"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[0,0,0],[2,2,2]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 1: sx is greater than tolerance. |sx - reference| = 2e0"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[0,0,0],[0,2,2]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 1: sy is greater than tolerance. |sy - reference| = 2e0"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2]],"stresses":[[[0,0,0],[0,0,2]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 1: sxy is greater than tolerance. |sxy - reference| = 2e0"
        );
        Ok(())
    }

    #[test]
    fn compare_stresses_captures_errors_3d() -> Result<(), StrError> {
        let two_dim = false;
        let mut equations = EquationNumbers::new(1);
        equations.activate_equation(0, Dof::Ux);
        equations.activate_equation(0, Dof::Uy);
        equations.activate_equation(0, Dof::Uz);
        let neq = equations.nequation();
        let sigma = StateStress {
            internal_values: Vec::new(),
            sigma: Tensor2::from_matrix(&[[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]], true, two_dim)?,
        };
        let mut state = State {
            elements: vec![StateElement {
                seepage: Vec::new(),
                stress: vec![sigma.clone(), sigma.clone()],
            }],
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
        };
        state.system_xx[0] = 3.0;
        state.system_xx[1] = 4.0;
        state.system_xx[2] = 6.0;

        let val = Validator::from_str(r#"{"steps":[{"disp":[[3,4,6]],"stresses":[[[0,0,0]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 0: reference stress has incompatible number of components. 3(wrong) != 4"
        );

        let val = Validator::from_str(r#"{"steps":[{"disp":[[3,4,6]],"stresses":[[[0,0,0,2]]]}]}"#)?;
        assert_eq!(
            val.compare_stresses(0, &state, two_dim),
            "element 0: integration point 0: sz is greater than tolerance. |sz - reference| = 2e0"
        );
        Ok(())
    }

    #[test]
    fn compare_displacements_and_stresses_work_with_no_data() -> Result<(), StrError> {
        let two_dim = true;
        let equations = EquationNumbers::new(0);
        let state = State::new_empty();
        let val = Validator::from_str(r#"{"steps":[{}]}"#)?;
        assert_eq!(val.compare_displacements(0, &state, &equations, two_dim), "OK");
        assert_eq!(val.compare_stresses(0, &state, two_dim), "OK");
        Ok(())
    }

    #[test]
    fn compare_displacements_and_stress_work_2d() -> Result<(), StrError> {
        /* Example 1.6 from [@bhatti] page 32

         Solid bracket with thickness = 0.25

                      1    load                connectivity:
         y=2.0  fixed *'-,__                    eid : vertices
                      |     '-,_  3   load        0 :  0, 2, 3
         y=1.5 - - -  |        ,'*-,__            1 :  3, 1, 0
                      |  1   ,'  |    '-,_  5     2 :  2, 4, 5
         y=1.0 - - -  |    ,'    |  3   ,-'*      3 :  5, 3, 2
                      |  ,'  0   |   ,-'   |
                      |,'        |,-'   2  |   constraints:
         y=0.0  fixed *----------*---------*     fixed on x and y
                      0          2         4
                     x=0.0     x=2.0     x=4.0

        # References

        [@bhatti] Bhatti, M.A. (2005) Fundamental Finite Element Analysis
                  and Applications, Wiley, 700p.
        */

        let mut equations = EquationNumbers::new(6);
        for point_id in 0..6 {
            equations.activate_equation(point_id, Dof::Ux);
            equations.activate_equation(point_id, Dof::Uy);
        }

        let neq = equations.nequation();
        let two_dim = true;

        let state_element_0 = StateElement {
            seepage: Vec::new(),
            stress: vec![StateStress {
                sigma: Tensor2::from_matrix(
                    &[
                        [-5.283090599362460e+01, -1.128984616188524e+01, 0.0],
                        [-1.128984616188524e+01, -5.272560566371797e+00, 0.0],
                        [0.0, 0.0, 0.0],
                    ],
                    true,
                    two_dim,
                )?,
                internal_values: Vec::new(),
            }],
        };

        let state_element_1 = StateElement {
            seepage: Vec::new(),
            stress: vec![StateStress {
                sigma: Tensor2::from_matrix(
                    &[
                        [2.462317949521848e+01, -5.153261537858599e+01, 0.0],
                        [-5.153261537858599e+01, 4.924635899043697e+00, 0.0],
                        [0.0, 0.0, 0.0],
                    ],
                    true,
                    two_dim,
                )?,
                internal_values: Vec::new(),
            }],
        };

        let state = State {
            elements: vec![state_element_0, state_element_1],
            system_xx: Vector::from(&[
                0.000000000000000e+00,
                0.000000000000000e+00,
                0.000000000000000e+00,
                0.000000000000000e+00,
                -1.035527877607004e-02,
                -2.552969847657423e-02,
                4.727650463081949e-03,
                -2.473565538172127e-02,
                -1.313941349422282e-02,
                -5.549310752960183e-02,
                8.389015766816341e-05,
                -5.556637423271112e-02,
            ]),
            system_yy: Vector::new(neq),
        };

        let mut val = Validator::read_json("./data/validation/bhatti_1_6.json")?;
        val.tol_displacement = 1e-14;
        val.tol_stress = 1e-14;
        assert_eq!(val.compare_displacements(0, &state, &equations, two_dim), "OK");
        assert_eq!(val.compare_stresses(0, &state, two_dim), "OK");
        Ok(())
    }

    #[test]
    fn compare_displacements_and_stress_work_3d() -> Result<(), StrError> {
        let mut equations = EquationNumbers::new(1);
        equations.activate_equation(0, Dof::Ux);
        equations.activate_equation(0, Dof::Uy);
        equations.activate_equation(0, Dof::Uz);

        let neq = equations.nequation();
        let two_dim = false;

        let state_element_0 = StateElement {
            seepage: Vec::new(),
            stress: vec![StateStress {
                sigma: Tensor2::from_matrix(&[[1.1, 1.2, 1.3], [1.2, 2.2, 2.3], [1.3, 2.3, 3.3]], true, two_dim)?,
                internal_values: Vec::new(),
            }],
        };

        let state = State {
            elements: vec![state_element_0],
            system_xx: Vector::from(&[1.0, 2.0, 3.0]),
            system_yy: Vector::new(neq),
        };

        let val = Validator::from_str(r#"{"steps":[{"disp":[[1,2,3]],"stresses":[[[1.1,2.2,1.2,3.3]]]}]}"#)?;
        assert_eq!(val.compare_displacements(0, &state, &equations, two_dim), "OK");
        assert_eq!(val.compare_stresses(0, &state, two_dim), "OK");
        Ok(())
    }
}
