use crate::base::{Config, Dof};
use crate::fem::{FemBase, FemState};
use crate::StrError;
use gemlab::mesh::{CellId, Mesh, PointId};
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a tensor-valued quantity at a Gauss point along time
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct TemporalTensor {
    pub txx: Vec<f64>,
    pub tyy: Vec<f64>,
    pub tzz: Vec<f64>,
    pub txy: Vec<f64>,
    pub tyz: Vec<f64>,
    pub tzx: Vec<f64>,
}

impl TemporalTensor {
    /// Allocates a new instance
    pub fn new() -> Self {
        TemporalTensor {
            txx: Vec::new(),
            tyy: Vec::new(),
            tzz: Vec::new(),
            txy: Vec::new(),
            tyz: Vec::new(),
            tzx: Vec::new(),
        }
    }
}

/// Assists in generating output files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemResults {
    /// Number of files written
    counter: usize,

    /// Indices of the output files
    pub indices: Vec<usize>,

    /// Real simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Step number for selected points and cells
    pub sel_step: Vec<usize>,

    /// Time for selected points and cells
    pub sel_time: Vec<f64>,

    /// Loading factors for selected points and cells
    pub sel_lambda: Vec<f64>,

    /// DOF values at selected points along time
    ///
    /// Maps "PointId,Dof" to an array with the values along time
    sel_dof: HashMap<String, Vec<f64>>,

    /// Stresses at selected integration points along time
    ///
    /// The results at the first integration point are saved only.
    sel_stress: HashMap<CellId, TemporalTensor>,

    /// Strains at selected integration points along time
    ///
    /// The results at the first integration point are saved only.
    sel_strain: HashMap<CellId, TemporalTensor>,
}

impl FemResults {
    /// Allocates a new instance with deactivated generation of files
    pub fn new(mesh: &Mesh, base: &FemBase, config: &Config) -> Result<Self, StrError> {
        if config.out_files {
            // create directory
            fs::create_dir_all(&config.out_dir).map_err(|_| "cannot create output directory")?;

            // write the mesh
            mesh.write(&format!("{}/{}-mesh.msh", config.out_dir, config.out_fn_stem))?;

            // write the FEM base
            base.write_json(&format!("{}/{}-base.json", config.out_dir, config.out_fn_stem))?;
        }
        Ok(FemResults {
            counter: 0,
            indices: Vec::new(),
            times: Vec::new(),
            sel_time: Vec::new(),
            sel_step: Vec::new(),
            sel_lambda: Vec::new(),
            sel_dof: HashMap::new(),
            sel_stress: HashMap::new(),
            sel_strain: HashMap::new(),
        })
    }

    /// Returns the temporal output of DOF values at selected points
    pub fn get_dof(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        let key = format!("{:?},{:?}", point_id, dof);
        self.sel_dof.get(&key)
    }

    /// Returns the temporal output of stresses at selected integration points
    pub fn get_stress(&self, cell_id: CellId) -> Option<&TemporalTensor> {
        self.sel_stress.get(&cell_id)
    }

    /// Returns the temporal output of strains at selected integration points
    pub fn get_strain(&self, cell_id: CellId) -> Option<&TemporalTensor> {
        self.sel_strain.get(&cell_id)
    }

    /// Reads a JSON file containing this struct
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open JSON file")?;
        let buffered = BufReader::new(data);
        let results = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
        Ok(results)
    }

    /// Writes a JSON file with this struct
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let mut file = File::create(&path).map_err(|_| "cannot create JSON file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write JSON file")?;
        Ok(())
    }

    /// Writes the current FEM state to a file
    pub(crate) fn write_state(&mut self, config: &Config, state: &FemState) -> Result<(), StrError> {
        if config.out_files {
            // save the state
            state.write_json(&format!(
                "{}/{}-{:0>20}.json",
                config.out_dir, config.out_fn_stem, self.counter
            ))?;

            // update counters
            self.indices.push(self.counter);
            self.times.push(state.time);
            self.counter += 1;
        }
        Ok(())
    }

    /// Writes this struct to a file
    pub(crate) fn write_self(&self, config: &Config) -> Result<(), StrError> {
        if config.out_files {
            self.write_json(&format!("{}/{}.json", config.out_dir, config.out_fn_stem))?;
        }
        Ok(())
    }

    /// Saves the results at selected nodes and integration points
    pub(crate) fn save_selected(&mut self, config: &Config, base: &FemBase, state: &FemState) -> Result<(), StrError> {
        if config.out_has_selected {
            // step, time, and lambda
            self.sel_step.push(state.step);
            self.sel_time.push(state.time);
            self.sel_lambda.push(state.lambda);

            // DOFs
            for (point_id, dof) in config.out_dof.iter() {
                if let Some(eq) = base.dofs.eq(*point_id, *dof).ok() {
                    let key = format!("{:?},{:?}", point_id, dof);
                    self.sel_dof.entry(key).or_insert(Vec::new()).push(state.u[eq]);
                }
            }

            // stresses
            for cell_id in config.out_stress.iter() {
                if let Some(s) = state.gauss[*cell_id].stress(0).ok() {
                    let r = self.sel_stress.entry(*cell_id).or_insert(TemporalTensor::new());
                    r.txx.push(s.get(0, 0));
                    r.tyy.push(s.get(1, 1));
                    r.tzz.push(s.get(2, 2));
                    r.txy.push(s.get(0, 1));
                    if s.dim() > 4 {
                        r.tyz.push(s.get(1, 2));
                        r.tzx.push(s.get(2, 0));
                    }
                }
            }

            // strains
            for cell_id in config.out_strain.iter() {
                if let Some(s) = state.gauss[*cell_id].strain(0).ok() {
                    let r = self.sel_strain.entry(*cell_id).or_insert(TemporalTensor::new());
                    r.txx.push(s.get(0, 0));
                    r.tyy.push(s.get(1, 1));
                    r.tzz.push(s.get(2, 2));
                    r.txy.push(s.get(0, 1));
                    if s.dim() > 4 {
                        r.tyz.push(s.get(1, 2));
                        r.tzx.push(s.get(2, 0));
                    }
                }
            }
        }
        Ok(())
    }
}
