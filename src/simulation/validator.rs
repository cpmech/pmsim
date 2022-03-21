use crate::StrError;
use serde::{Deserialize, Serialize};
use serde_json;

/// Holds results from iterations
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ValidatorIteration {
    /// iteration number
    #[serde(rename(deserialize = "it"))]
    pub iteration: usize,

    /// relative residual
    #[serde(rename(deserialize = "resrel"))]
    pub relative_residual: f64,

    /// absolute residual
    #[serde(rename(deserialize = "resid"))]
    pub absolute_residual: f64,
}

/// Holds numerical results
///
/// Stresses: examples:
///   2D solid: [sx, sy, sxy, sz]
///   3D beam:  [M22, M11, Mtt, Vs, Vr]
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ValidatorResults {
    /// all stiffness matrices (nele,nu,nu)
    #[serde(default)]
    #[serde(rename(deserialize = "Kmats"))]
    pub kk_matrices: Vec<Vec<Vec<f64>>>,

    /// displacements at nodes (npoint,ndim)
    #[serde(default)]
    #[serde(rename(deserialize = "disp"))]
    pub displacements: Vec<Vec<f64>>,

    /// all stresses @ all ips (nele,nip,nsigma)
    #[serde(default)]
    pub stresses: Vec<Vec<Vec<f64>>>,

    /// load factor
    #[serde(default)]
    #[serde(rename(deserialize = "loadfactor"))]
    pub load_factor: f64,

    /// iterations data
    #[serde(default)]
    pub iterations: Vec<ValidatorIteration>,
}

/// Holds results for comparisons (checking/tests)
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct Validator {
    /// results for each time-step output
    pub results: Vec<ValidatorResults>,
}

impl Validator {
    /// Returns a new Validator from JSON string
    pub fn from_json(cmp: &str) -> Result<Self, StrError> {
        let res: Validator = serde_json::from_str(&cmp).map_err(|op| {
            println!("ERROR: {}", op);
            return "serde_json failed";
        })?;
        Ok(res)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Validator;
    use crate::StrError;
    use russell_chk::assert_vec_approx_eq;

    #[test]
    fn clone_and_serialize_work() -> Result<(), StrError> {
        let val = Validator::from_json(
            "{ \"results\":\n\
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
        assert_eq!(val.results.len(), 1);
        let res = &val.results[0];
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
        assert_eq!(cloned.results.len(), 1);
        let c_res = &cloned.results[0];
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
            "{ \"results\":\n\
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
        assert_eq!(val.results.len(), 3);
        assert_eq!(val.results[0].iterations.len(), 1);
        assert_eq!(val.results[1].iterations.len(), 2);
        assert_eq!(val.results[2].iterations.len(), 3);
        Ok(())
    }
}
