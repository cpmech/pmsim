use crate::base::{Dof, Equations, DEFAULT_OUT_DIR};
use crate::fem::{FemMesh, FemState};
use crate::StrError;
use gemlab::mesh::{At, Features, PointId};
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a summary of the generated files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemOutputSummary {
    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,

    /// Holds equation numbers (DOF numbers)
    pub equations: Option<Equations>,
}

/// Assists in the post-processing of results
pub struct FemOutput<'a> {
    /// Holds the FEM mesh, attributes, and DOF numbers
    fem: &'a FemMesh<'a>,

    /// Defines the filename stem
    filename_stem: Option<String>,

    /// Defines the output directory
    output_directory: String,

    /// Holds the count of files written
    output_count: usize,

    /// Holds the summary
    summary: FemOutputSummary,
}

impl FemOutputSummary {
    /// Reads a JSON file containing the summary
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let fem = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(fem);
        let summary = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
        Ok(summary)
    }

    /// Writes a JSON file with the summary
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }
}

impl<'a> FemOutput<'a> {
    /// Allocates new instance
    ///
    /// # Input
    ///
    /// * `fem` -- the FEM mesh, attributes, and DOF numbers
    /// * `filename_steam` -- the last part of the filename without extension, e.g., "my_simulation".
    ///   None means that no files will be written.
    /// * `output_directory` -- the directory to save the output files.
    ///   None means that the default directory will be used; see [DEFAULT_OUT_DIR]
    pub fn new(
        fem: &'a FemMesh,
        filename_stem: Option<String>,
        output_directory: Option<&str>,
    ) -> Result<Self, StrError> {
        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUT_DIR,
        };

        // create directory only if filename_stem is provided
        if let Some(_) = filename_stem {
            fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;
        }

        // summary
        let summary = match &filename_stem {
            Some(_) => FemOutputSummary {
                indices: Vec::new(),
                times: Vec::new(),
                equations: Some(fem.equations.clone()),
            },
            None => FemOutputSummary {
                indices: Vec::new(),
                times: Vec::new(),
                equations: None,
            },
        };

        // output
        Ok(FemOutput {
            fem,
            filename_stem,
            output_directory: out_dir.to_string(),
            output_count: 0,
            // callback,
            summary,
        })
    }

    /// Generates the filename path for the mesh file
    pub fn path_mesh(out_dir: &str, fn_stem: &str) -> String {
        format!("{}/{}-mesh.json", out_dir, fn_stem)
    }

    /// Generates the filename path for the summary file
    pub fn path_summary(out_dir: &str, fn_stem: &str) -> String {
        format!("{}/{}-summary.json", out_dir, fn_stem)
    }

    /// Generates the filename path for the state files
    pub fn path_state(out_dir: &str, fn_stem: &str, index: usize) -> String {
        format!("{}/{}-{:0>20}.json", out_dir, fn_stem, index)
    }

    /// Writes the current FEM state to a file
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write(&mut self, state: &mut FemState) -> Result<(), StrError> {
        if let Some(fn_stem) = &self.filename_stem {
            // save the mesh
            if self.output_count == 0 {
                let path = &FemOutput::path_mesh(&self.output_directory, fn_stem);
                self.fem.mesh.write_json(&path)?;
            }

            // save the state
            let path = FemOutput::path_state(&self.output_directory, fn_stem, self.output_count);
            state.write_json(&path)?;

            // update summary
            self.summary.indices.push(self.output_count);
            self.summary.times.push(state.t);
            self.output_count += 1;
        }
        Ok(())
    }

    /// Writes the summary of generated files (at the end of the simulation)
    pub(crate) fn write_summary(&self) -> Result<(), StrError> {
        if let Some(fn_stem) = &self.filename_stem {
            let path = FemOutput::path_summary(&self.output_directory, fn_stem);
            self.summary.write_json(&path)?;
        }
        Ok(())
    }

    /// Extracts primary values along x
    ///
    /// Returns the DOF values (e.g., T) corresponding to the points with constant y coordinate
    ///
    /// **Important:** If you need values at points on the interior of the mesh,
    /// then you have to pass the Extract::All option when allocating a new Find instance.
    ///
    /// # Input
    ///
    /// * `state` -- state for which the {U} values are extracted
    /// * `dof` -- the desired DOF, e.g., T
    /// * `y` -- the constant elevation
    /// * `filter` -- fn(x) -> bool that returns true to **keep** the coordinate just found
    ///   (yields only the elements for which the closure returns true)
    ///
    /// # Output
    ///
    /// Returns `(ids, xx, dd)`, where:
    ///
    /// * `ids` -- contains the IDs of the points along x
    /// * `xx` -- are the x-coordinates
    /// * `dd` -- are the DOF values (e.g., temperature) along x and corresponding to the `ids` and `xx`
    pub fn values_along_x<F>(
        &self,
        features: &Features,
        state: &FemState,
        dof: Dof,
        y: f64,
        filter: F,
    ) -> Result<(Vec<PointId>, Vec<f64>, Vec<f64>), StrError>
    where
        F: FnMut(&[f64]) -> bool,
    {
        // find points and sort by x-coordinates
        let point_ids = features.search_point_ids(At::Y(y), filter)?;
        let mut id_x_pairs: Vec<_> = point_ids
            .iter()
            .map(|id| (*id, self.fem.mesh.points[*id].coords[0]))
            .collect();
        id_x_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // extract dof values
        let dd: Vec<_> = id_x_pairs
            .iter()
            .map(|(id, _)| state.uu[self.fem.equations.eq(*id, dof).unwrap()])
            .collect();

        // unzip id_x_pairs
        let (ids, xx): (Vec<_>, Vec<_>) = id_x_pairs.iter().cloned().unzip();

        // results
        Ok((ids, xx, dd))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemOutput;
    use crate::base::{Config, Dof, Elem, ParamDiffusion};
    use crate::fem::{FemMesh, FemState};
    use gemlab::mesh::{Features, Samples};
    use gemlab::util::any_x;

    #[test]
    fn values_along_x_works() {
        let mesh = Samples::one_tri6();
        let features = Features::new(&mesh, false);
        let p1 = ParamDiffusion::sample();
        let fem = FemMesh::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&fem, &config).unwrap();
        state.uu[0] = 1.0;
        state.uu[1] = 2.0;
        state.uu[2] = 3.0;
        state.uu[3] = 4.0;
        state.uu[4] = 5.0;
        state.uu[5] = 6.0;
        let output = FemOutput::new(&fem, None, None).unwrap();
        let (ids, xx, dd) = output.values_along_x(&features, &state, Dof::T, 0.0, any_x).unwrap();
        assert_eq!(ids, &[0, 3, 1]);
        assert_eq!(xx, &[0.0, 0.5, 1.0]);
        assert_eq!(dd, &[1.0, 4.0, 2.0]);
    }
}
