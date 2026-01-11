use crate::base::{Config, Dof, Schema};
use crate::fem::FemState;
use crate::material::LocalState;
use crate::StrError;
use gemlab::mesh::{CellId, Mesh, PointId};
use russell_lab::Vector;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Assists in generating output files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub(crate) struct OutputFiles {
    /// Number of files written
    counter: usize,

    /// Indices of the output files
    indices: Vec<usize>,

    /// Real simulation times corresponding to each output file
    times: Vec<f64>,

    /// Time for selected points and cells
    sel_time: Vec<f64>,

    /// Loading factors for selected points and cells
    sel_lambda: Vec<f64>,

    /// DOF values at selected points along time
    ///
    /// Maps "PointId,Dof" to an array with the values along time
    sel_dof: HashMap<String, Vec<f64>>,

    /// Flux vectors at selected integration points along time
    ///
    /// Note: Only the results at the first integration point are saved.
    sel_local_flux: HashMap<CellId, Vec<Vector>>,

    /// LocalState at selected integration points along time
    ///
    /// Note: Only the results at the first integration point are saved.
    sel_local_state: HashMap<CellId, Vec<LocalState>>,
}

impl OutputFiles {
    /// Allocates a new instance with deactivated generation of files
    pub fn new(mesh: &Mesh, schema: &Schema, config: &Config) -> Result<Self, StrError> {
        if config.out_files {
            // create directory
            fs::create_dir_all(&config.out_dir).map_err(|_| "cannot create output directory")?;

            // write the mesh
            mesh.write(&format!("{}/{}-mesh.msh", config.out_dir, config.out_fn_stem))?;

            // write the FEM base
            schema.write_json(&format!("{}/{}-schema.json", config.out_dir, config.out_fn_stem))?;
        }
        Ok(OutputFiles {
            counter: 0,
            indices: Vec::new(),
            times: Vec::new(),
            sel_time: Vec::new(),
            sel_lambda: Vec::new(),
            sel_dof: HashMap::new(),
            sel_local_flux: HashMap::new(),
            sel_local_state: HashMap::new(),
        })
    }

    /// Returns the number of files written
    pub fn n_files(&self) -> usize {
        self.counter
    }

    /// Returns the indices of the output files
    pub fn get_indices(&self) -> &Vec<usize> {
        &self.indices
    }

    /// Returns the real simulation times corresponding to each output file
    pub fn get_times(&self) -> &Vec<f64> {
        &self.times
    }

    /// Returns the temporal output of DOF values at selected points
    pub fn get_selected_dof(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        let key = format!("{:?},{:?}", point_id, dof);
        self.sel_dof.get(&key)
    }

    /// Returns the temporal output of flux vectors at the first integration point of selected cells
    pub fn get_selected_local_fluxes(&self, cell_id: CellId) -> Option<&Vec<Vector>> {
        self.sel_local_flux.get(&cell_id)
    }

    /// Returns the temporal output of stresses at the first integration point of selected cells
    pub fn get_selected_local_state(&self, cell_id: CellId) -> Option<&Vec<LocalState>> {
        self.sel_local_state.get(&cell_id)
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

    /// Executes the output
    pub(crate) fn execute(&mut self, schema: &Schema, config: &Config, state: &FemState) -> Result<(), StrError> {
        if config.out_files {
            // save the state
            state.write_json(&format!(
                "{}/{}-{}.json",
                config.out_dir, config.out_fn_stem, self.counter
            ))?;

            // update counters
            self.indices.push(self.counter);
            self.times.push(state.time);
            self.counter += 1;
        }
        if config.out_has_selected {
            // step, time, and lambda
            self.sel_time.push(state.time);
            self.sel_lambda.push(state.lambda);

            // DOFs
            for (point_id, dof) in config.out_dof.iter() {
                if schema.has_dof(*point_id, *dof)? {
                    let eq = schema.get_eq(*point_id, *dof)?;
                    let key = format!("{:?},{:?}", point_id, dof);
                    self.sel_dof.entry(key).or_insert(Vec::new()).push(state.u[eq]);
                }
            }

            // local states
            for cell_id in config.out_local_state.iter() {
                if let Some(w) = state.gauss[*cell_id].get_flux_vector(0).ok() {
                    self.sel_local_flux
                        .entry(*cell_id)
                        .or_insert(Vec::new())
                        .push(w.clone());
                }
                if let Some(s) = state.gauss[*cell_id].get_local_state(0).ok() {
                    self.sel_local_state
                        .entry(*cell_id)
                        .or_insert(Vec::new())
                        .push(s.clone());
                }
            }
        }
        Ok(())
    }

    /// Stops the output
    pub(crate) fn stop(&self, config: &Config) -> Result<(), StrError> {
        if config.out_files {
            self.write_json(&format!("{}/{}.json", config.out_dir, config.out_fn_stem))?;
        }
        Ok(())
    }
}
