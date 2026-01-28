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

    /// Total number of equations
    neq_total: usize,

    /// Number of prescribed equations
    neq_presc: usize,

    /// Time for selected points and cells
    sel_time: Vec<f64>,

    /// Loading factors for selected points and cells
    sel_lambda: Vec<f64>,

    /// U components at selected points along time
    ///
    /// Maps `equation` to an array with the values along time
    sel_uu_comp: HashMap<usize, Vec<f64>>,

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
    pub fn new(mesh: &Mesh, schema: &Schema, config: &Config, neq_presc: usize) -> Result<Self, StrError> {
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
            neq_total: schema.get_neq()?,
            neq_presc,
            sel_time: Vec::new(),
            sel_lambda: Vec::new(),
            sel_uu_comp: HashMap::new(),
            sel_local_flux: HashMap::new(),
            sel_local_state: HashMap::new(),
        })
    }

    /// Starts the output
    pub fn start(&mut self) {
        self.counter = 0;
        self.indices.clear();
        self.times.clear();
        self.sel_time.clear();
        self.sel_lambda.clear();
        self.sel_uu_comp.clear();
        self.sel_local_flux.clear();
        self.sel_local_state.clear();
    }

    /// Returns the number of files written
    pub fn n_files(&self) -> usize {
        self.counter
    }

    /// Returns the total number of equations
    pub fn neq_total(&self) -> usize {
        self.neq_total
    }

    /// Returns the number of prescribed equations
    pub fn neq_presc(&self) -> usize {
        self.neq_presc
    }

    /// Returns the indices of the output files
    pub fn get_indices(&self) -> &Vec<usize> {
        &self.indices
    }

    /// Returns the real simulation times corresponding to each output file
    pub fn get_times(&self) -> &Vec<f64> {
        &self.times
    }

    /// Returns the temporal output of a selected U component
    pub fn get_selected_uu_comp(&self, point_id: PointId, dof: Dof, schema: &Schema) -> Option<&Vec<f64>> {
        match schema.get_eq(point_id, dof) {
            Ok(eq) => self.sel_uu_comp.get(&eq),
            Err(_) => None,
        }
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

            // U components
            for (point_id, dof) in config.out_uu_comp.iter() {
                if schema.has_dof(*point_id, *dof)? {
                    let eq = schema.get_eq(*point_id, *dof)?;
                    self.sel_uu_comp.entry(eq).or_insert(Vec::new()).push(state.u[eq]);
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
