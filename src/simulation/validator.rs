use super::{Dof, EquationNumbers, State};
use crate::StrError;
use russell_tensor::SQRT_2;
use serde::{Deserialize, Serialize};
use serde_json;
use std::ffi::OsStr;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

/// Holds results from iterations
#[derive(Clone, Debug, Deserialize, Serialize)]
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
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ValidatorResults {
    /// All stiffness matrices (nele,nu,nu)
    #[serde(default)]
    #[serde(rename(deserialize = "Kmats"))]
    pub kk_matrices: Vec<Vec<Vec<f64>>>,

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
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Validator {
    /// Results from output-, time-, or load- steps
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
    pub fn from_json(cmp: &str) -> Result<Self, StrError> {
        let res = serde_json::from_str(&cmp).map_err(|op| {
            println!("ERROR: {}", op);
            return "serde_json failed";
        })?;
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
        let res = serde_json::from_reader(reader).map_err(|op| {
            println!("ERROR: {}", op);
            return "serde_json failed";
        })?;
        Ok(res)
    }

    /// Compares state against results at a fixed time- or load- step
    pub fn compare_state(&self, step: usize, state: &State, equations: &EquationNumbers, two_dim: bool) -> String {
        if step >= self.steps.len() {
            return "results for the step are not available".to_string();
        }
        let cmp = &self.steps[step];

        // displacements
        let npoint = equations.npoint();
        for point_id in 0..npoint {
            if point_id >= cmp.displacements.len() {
                return format!("displacement of point #{} is not available", point_id);
            }
            let reference = &cmp.displacements[point_id];
            if reference.len() < 2 {
                return format!("reference (ux,uy) for point #{} is missing", point_id);
            }
            let ux = match equations.number(point_id, Dof::Ux) {
                Some(eq) => state.system_xx[eq],
                None => return format!("state does not have ux displacement of point #{}", point_id),
            };
            let uy = match equations.number(point_id, Dof::Uy) {
                Some(eq) => state.system_xx[eq],
                None => return format!("state does not have uy displacement of point #{}", point_id),
            };
            let diff_ux = f64::abs(ux - reference[0]);
            if diff_ux > self.tol_displacement {
                return format!(
                    "ux displacement of point #{} is greater than tolerance. |ux - reference| = {:e}",
                    point_id, diff_ux
                );
            }
            let diff_uy = f64::abs(uy - reference[1]);
            if diff_uy > self.tol_displacement {
                return format!(
                    "uy displacement of point #{} is greater than tolerance. |uy - reference| = {:e}",
                    point_id, diff_uy
                );
            }
            if !two_dim {
                if reference.len() != 3 {
                    return format!("reference uz for point #{} is missing", point_id);
                }
                let uz = match equations.number(point_id, Dof::Uz) {
                    Some(eq) => state.system_xx[eq],
                    None => return format!("state does not have uz displacement of point #{}", point_id),
                };
                let diff_uz = f64::abs(uz - reference[2]);
                if diff_uz > self.tol_displacement {
                    return format!(
                        "uz displacement of point #{} is greater than tolerance. |uz - reference| = {:e}",
                        point_id, diff_uz
                    );
                }
            }
        }

        // stresses
        let nele = state.elements.len();
        for element_id in 0..nele {
            let element = &state.elements[element_id];
            let n_integ_point = element.stress.len();
            if n_integ_point > 0 {
                if element_id >= cmp.stresses.len() {
                    return format!("stresses for element #{} are not available", element_id);
                }
                for index_ip in 0..n_integ_point {
                    if index_ip >= cmp.stresses[element_id].len() {
                        return format!(
                            "stress at integration point #{} of element #{} is missing",
                            index_ip, element_id
                        );
                    }
                    let reference = &cmp.stresses[element_id][index_ip];
                    if reference.len() < 3 {
                        return format!(
                            "reference (sx,sy,sxy) for element #{} at ip #{} is missing",
                            element_id, index_ip
                        );
                    }
                    let sigma = &element.stress[index_ip].stress;
                    let sx = sigma.vec[0];
                    let sy = sigma.vec[1];
                    let sxy = sigma.vec[3] / SQRT_2;
                    let diff_sx = f64::abs(sx - reference[0]);
                    if diff_sx > self.tol_stress {
                        return format!(
                            "sx stress of element #{} at ip #{} is greater than tolerance. |sx - reference| = {:e}",
                            element_id, index_ip, diff_sx
                        );
                    }
                    let diff_sy = f64::abs(sy - reference[1]);
                    if diff_sy > self.tol_stress {
                        return format!(
                            "sy stress of element #{} at ip #{} is greater than tolerance. |sy - reference| = {:e}",
                            element_id, index_ip, diff_sy
                        );
                    }
                    let diff_sxy = f64::abs(sxy - reference[2]);
                    if diff_sxy > self.tol_stress {
                        return format!(
                            "sxy stress of element #{} at ip #{} is greater than tolerance. |sxy - reference| = {:e}",
                            element_id, index_ip, diff_sxy
                        );
                    }
                    if !two_dim {
                        if reference.len() != 4 {
                            return format!(
                                "reference sz for element #{} at ip #{} is missing",
                                element_id, index_ip
                            );
                        }
                        let sz = sigma.vec[2];
                        let diff_sz = f64::abs(sz - reference[3]);
                        if diff_sz > self.tol_stress {
                            return format!(
                                "sy stress of element #{} at ip #{} is greater than tolerance. |sz - reference| = {:e}",
                                element_id, index_ip, diff_sz
                            );
                        }
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
    use super::Validator;
    use crate::simulation::{Dof, EquationNumbers, State, StateElement, StateStress};
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Vector;
    use russell_tensor::Tensor2;

    #[test]
    fn clone_and_serialize_work() -> Result<(), StrError> {
        let val = Validator::from_json(
            "{ \"steps\":\n\
          [\n\
            {\n\
              \"Kmats\": [\n\
                [\n\
                  [1, 2],\n\
                  [3, 4]\n\
                ],\n\
                [\n\
                  [10, 20],\n\
                  [30, 40]\n\
                ]\n\
              ],\n\
              \"disp\": [\n\
                [11, 21],\n\
                [12, 22]\n\
              ],\n\
              \"stresses\": [\n\
                [\n\
                  [100, 101, 102, 103]\n\
                ],\n\
                [\n\
                  [200, 201, 202, 203]\n\
                ]\n\
              ]\n\
            }\n\
          ]\n\
        }",
        )?;
        assert_eq!(val.steps.len(), 1);
        let res = &val.steps[0];
        assert_eq!(res.kk_matrices.len(), 2);
        assert_vec_approx_eq!(res.kk_matrices[0][0], [1.0, 2.0], 1e-15);
        assert_vec_approx_eq!(res.kk_matrices[0][1], [3.0, 4.0], 1e-15);
        assert_vec_approx_eq!(res.kk_matrices[1][0], [10.0, 20.0], 1e-15);
        assert_vec_approx_eq!(res.kk_matrices[1][1], [30.0, 40.0], 1e-15);
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
        assert_vec_approx_eq!(c_res.kk_matrices[0][0], [1.0, 2.0], 1e-15);
        assert_vec_approx_eq!(c_res.kk_matrices[0][1], [3.0, 4.0], 1e-15);
        assert_vec_approx_eq!(c_res.kk_matrices[1][0], [10.0, 20.0], 1e-15);
        assert_vec_approx_eq!(c_res.kk_matrices[1][1], [30.0, 40.0], 1e-15);
        assert_eq!(c_res.displacements.len(), 2);
        assert_vec_approx_eq!(c_res.displacements[0], [11.0, 21.0], 1e-15);
        assert_vec_approx_eq!(c_res.displacements[1], [12.0, 22.0], 1e-15);
        assert_eq!(c_res.stresses.len(), 2);
        assert_vec_approx_eq!(c_res.stresses[0][0], [100.0, 101.0, 102.0, 103.0], 1e-15);
        assert_vec_approx_eq!(c_res.stresses[1][0], [200.0, 201.0, 202.0, 203.0], 1e-15);

        assert_eq!(format!("{:?}", res),"ValidatorResults { kk_matrices: [[[1.0, 2.0], [3.0, 4.0]], [[10.0, 20.0], [30.0, 40.0]]], displacements: [[11.0, 21.0], [12.0, 22.0]], stresses: [[[100.0, 101.0, 102.0, 103.0]], [[200.0, 201.0, 202.0, 203.0]]], load_factor: 0.0, iterations: [] }");

        Ok(())
    }

    #[test]
    fn from_json_works() -> Result<(), StrError> {
        let val = Validator::from_json(
            "{ \"steps\":\n\
          [\n\
            {\n\
              \"loadfactor\" : 1,\n\
              \"iterations\" : [\n\
                { \"it\":1, \"resrel\": 1, \"resid\": 2}\n\
              ],\n\
              \"disp\" : [\n\
                [1, 2],\n\
                [1, 2]\n\
              ],\n\
              \"stresses\" : [\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ],\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ]\n\
              ]\n\
            },\n\
            {\n\
              \"loadfactor\" : 2,\n\
              \"iterations\" : [\n\
                { \"it\":1, \"resrel\": 1, \"resid\": 2},\n\
                { \"it\":2, \"resrel\": 1, \"resid\": 2}\n\
              ],\n\
              \"disp\" : [\n\
                [1, 2],\n\
                [1, 2]\n\
              ],\n\
              \"stresses\" : [\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ],\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ]\n\
              ]\n\
            },\n\
            {\n\
              \"loadfactor\" : 3,\n\
              \"iterations\" : [\n\
                { \"it\":1, \"resrel\": 1, \"resid\": 2},\n\
                { \"it\":2, \"resrel\": 1, \"resid\": 2},\n\
                { \"it\":3, \"resrel\": 1, \"resid\": 2}\n\
              ],\n\
              \"disp\" : [\n\
                [1, 2],\n\
                [1, 2]\n\
              ],\n\
              \"stresses\" : [\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ],\n\
                [\n\
                  [1, 2, 3, 4],\n\
                  [1, 2, 3, 4]\n\
                ]\n\
              ]\n\
            }\n\
          ]\n\
        }",
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
    fn compare_state_works() -> Result<(), StrError> {
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
                stress: Tensor2::from_matrix(
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

        let state = State {
            elements: vec![state_element_0],
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

        let res = val.compare_state(0, &state, &equations, two_dim);
        assert_eq!(res, "OK");

        Ok(())
    }
}
