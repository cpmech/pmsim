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

    /// Real simulation times (time) or loading increments (lambda)
    ///
    /// Time is used for transient/dynamic analyses, while lambda is used for steady/static analyses
    stations: Vec<f64>,

    /// Total number of DOFs
    ndof: usize,

    /// Number of prescribed DOFs
    np: usize,

    /// History (time or lambda) of U components at selected points
    ///
    /// Maps `equation` to an array with the values along station
    history_uu_comp: HashMap<usize, Vec<f64>>,

    /// History (time or lambda) of Y (internal forces) components at selected points
    ///
    /// Maps `equation` to an array with the values along station
    history_yy_comp: HashMap<usize, Vec<f64>>,

    /// History (time or lambda) of flux vectors at selected integration points
    ///
    /// Note: Only the results at the first integration point are saved.
    history_local_flux: HashMap<CellId, Vec<Vector>>,

    /// History (time or lambda) of LocalState at selected integration points
    ///
    /// Note: Only the results at the first integration point are saved.
    history_local_state: HashMap<CellId, Vec<LocalState>>,
}

impl OutputFiles {
    /// Allocates a new instance with deactivated generation of files
    pub fn new(mesh: &Mesh, schema: &Schema, config: &Config, np: usize) -> Result<Self, StrError> {
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
            stations: Vec::new(),
            ndof: schema.ndof()?,
            np,
            history_uu_comp: HashMap::new(),
            history_yy_comp: HashMap::new(),
            history_local_flux: HashMap::new(),
            history_local_state: HashMap::new(),
        })
    }

    /// Returns the number of files written
    pub fn nfile(&self) -> usize {
        self.counter
    }

    /// Returns the total number of degrees of freedom (DOF)
    pub fn ndof(&self) -> usize {
        self.ndof
    }

    /// Returns the number of prescribed degrees of freedom (DOF)
    pub fn np(&self) -> usize {
        self.np
    }

    /// Returns the real simulation times (time) or loading increments (lambda)
    ///
    /// Time is used for transient/dynamic analyses, while lambda is used for steady/static analyses
    pub fn stations(&self) -> &Vec<f64> {
        &self.stations
    }

    /// Returns the history (time or lambda) of U components at selected points
    pub fn history_uu_comp(&self, point_id: PointId, dof: Dof, schema: &Schema) -> Option<&Vec<f64>> {
        match schema.dof_number(point_id, dof) {
            Ok(i) => self.history_uu_comp.get(&i),
            Err(_) => None,
        }
    }

    /// Returns the history (time or lambda) of Y (internal forces) components at selected points
    pub fn history_yy_comp(&self, point_id: PointId, dof: Dof, schema: &Schema) -> Option<&Vec<f64>> {
        match schema.dof_number(point_id, dof) {
            Ok(i) => self.history_yy_comp.get(&i),
            Err(_) => None,
        }
    }

    /// Returns the history (time or lambda) of flux vectors at selected integration points
    pub fn history_local_flux(&self, cell_id: CellId) -> Option<&Vec<Vector>> {
        self.history_local_flux.get(&cell_id)
    }

    /// Returns the history (time or lambda) of LocalState at selected integration points
    pub fn history_local_state(&self, cell_id: CellId) -> Option<&Vec<LocalState>> {
        self.history_local_state.get(&cell_id)
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
    pub(crate) fn execute(
        &mut self,
        schema: &Schema,
        config: &Config,
        state: &FemState,
        yy: &Vector,
    ) -> Result<(), StrError> {
        if config.out_files || config.out_history {
            if config.transient || config.dynamics {
                self.stations.push(state.time);
            } else {
                self.stations.push(state.lambda);
            }
        }
        if config.out_files {
            // save the state
            state.write_json(&format!(
                "{}/{}-{}.json",
                config.out_dir, config.out_fn_stem, self.counter
            ))?;
            self.counter += 1;
        }
        if config.out_history {
            // U components
            for (point_id, dof) in config.out_history_uu_comp.iter() {
                if schema.has_dof(*point_id, *dof)? {
                    let i = schema.dof_number(*point_id, *dof)?;
                    self.history_uu_comp.entry(i).or_insert(Vec::new()).push(state.uu[i]);
                }
            }

            // Y components
            for (point_id, dof) in config.out_history_yy_comp.iter() {
                if schema.has_dof(*point_id, *dof)? {
                    let i = schema.dof_number(*point_id, *dof)?;
                    self.history_yy_comp.entry(i).or_insert(Vec::new()).push(yy[i]);
                }
            }

            // local fluxes
            for cell_id in config.out_history_local_flux.iter() {
                if let Some(w) = state.gauss[*cell_id].get_flux_vector(0).ok() {
                    self.history_local_flux
                        .entry(*cell_id)
                        .or_insert(Vec::new())
                        .push(w.clone());
                }
            }

            // local states
            for cell_id in config.out_history_local_state.iter() {
                if let Some(s) = state.gauss[*cell_id].get_local_state(0).ok() {
                    self.history_local_state
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
