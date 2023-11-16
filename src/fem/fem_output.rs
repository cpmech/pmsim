use crate::base::{Dof, DEFAULT_OUT_DIR};
use crate::fem::{Elements, FemInput, FemState};
use crate::StrError;
use gemlab::mesh::{At, Features, PointId};
use serde::{Deserialize, Serialize};
use serde_json;
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds a summary of the generated files
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemOutputSummary {
    /// Last part of the filename without extension
    pub filename_stem: String,

    /// Space dimension
    pub ndim: usize,

    /// Number of points in the mesh
    pub npoint: usize,

    /// Number of cells in the mesh
    pub ncell: usize,

    /// Holds a subset of the currently DOFs in use
    pub dofs_in_use: Vec<Dof>,

    /// Holds a subset of the currently DOFs in use that aren't displacement DOF
    pub not_displacement_dof: Vec<Dof>,

    /// Holds the indices of the output files
    pub indices: Vec<usize>,

    /// Holds the simulation times corresponding to each output file
    pub times: Vec<f64>,
}

/// Assists in the post-processing of results
pub struct FemOutput<'a> {
    input: &'a FemInput<'a>,
    filename_stem: Option<String>,
    output_directory: String,
    output_count: usize,
    callback: Option<fn(&FemState, usize)>,
    summary: FemOutputSummary,
}

impl FemOutputSummary {
    /// Allocates a new instance
    pub(crate) fn new(input: &FemInput, filename_stem: String) -> Self {
        // DOFs in use
        let mut enabled_dofs = HashSet::new();
        for map in &input.equations.all {
            for dof in map.keys() {
                enabled_dofs.insert(*dof);
            }
        }
        let mut dofs_in_use: Vec<_> = enabled_dofs.iter().copied().collect();
        dofs_in_use.sort();

        // DOFs that aren't displacement DOFs
        let not_displacement_dof: Vec<_> = dofs_in_use
            .iter()
            .filter(|&&dof| !(dof == Dof::Ux || dof == Dof::Uy || dof == Dof::Uz))
            .copied()
            .collect();

        // summary
        FemOutputSummary {
            filename_stem,
            ndim: input.mesh.ndim,
            npoint: input.mesh.points.len(),
            ncell: input.mesh.cells.len(),
            dofs_in_use,
            not_displacement_dof,
            indices: Vec::new(),
            times: Vec::new(),
        }
    }

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
        let input = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(input);
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
    /// * `input` -- the FEM input data
    /// * `filename_steam` -- the last part of the filename without extension, e.g., "my_simulation".
    ///   None means that no files will be written.
    /// * `output_directory` -- the directory to save the output files.
    ///   None means that the default directory will be used; see [DEFAULT_OUTPUT_DIR]
    /// * `callback` -- is a function to be called at each time-output.
    ///   Example use `Some(|state, count| { ... })` to perform some processing on state at time `state.t`.
    ///   `count` is the corresponding `output_count` used to generate the output file.
    pub fn new(
        input: &'a FemInput,
        filename_stem: Option<String>,
        output_directory: Option<&str>,
        callback: Option<fn(&FemState, usize)>,
    ) -> Result<Self, StrError> {
        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUT_DIR,
        };
        fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;

        // make a copy of the filename_stem
        let fn_stem = if let Some(f) = &filename_stem {
            f.clone()
        } else {
            "".to_string()
        };

        // output
        Ok(FemOutput {
            input,
            filename_stem,
            output_directory: out_dir.to_string(),
            output_count: 0,
            callback,
            summary: FemOutputSummary::new(input, fn_stem),
        })
    }

    /// Writes the current FEM state to a binary file
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write(&mut self, state: &mut FemState, elements: &mut Elements) -> Result<(), StrError> {
        if let Some(fn_stem) = &self.filename_stem {
            elements.output_internal_values(state)?;
            let path = format!("{}/{}-{:0>20}.state", self.output_directory, fn_stem, self.output_count);
            state.write(&path)?;
            if let Some(callback) = self.callback {
                (callback)(state, self.output_count);
            }
            self.summary.indices.push(self.output_count);
            self.summary.times.push(state.t);
            self.output_count += 1;
        }
        Ok(())
    }

    /// Writes the summary of generated files (at the end of the simulation)
    pub(crate) fn write_summary(&self) -> Result<(), StrError> {
        if let Some(fn_stem) = &self.filename_stem {
            let path = format!("{}/{}-summary.json", self.output_directory, fn_stem);
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
        feat: &Features,
        state: &FemState,
        dof: Dof,
        y: f64,
        filter: F,
    ) -> Result<(Vec<PointId>, Vec<f64>, Vec<f64>), StrError>
    where
        F: FnMut(&[f64]) -> bool,
    {
        // find points and sort by x-coordinates
        let point_ids = feat.search_point_ids(At::Y(y), filter)?;
        let mut id_x_pairs: Vec<_> = point_ids
            .iter()
            .map(|id| (*id, self.input.mesh.points[*id].coords[0]))
            .collect();
        id_x_pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

        // extract dof values
        let dd: Vec<_> = id_x_pairs
            .iter()
            .map(|(id, _)| state.uu[self.input.equations.eq(*id, dof).unwrap()])
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
    use crate::base::{Config, Dof, Element, SampleParams};
    use crate::fem::{FemInput, FemState};
    use gemlab::mesh::{Features, Samples};
    use gemlab::util::any_x;

    #[test]
    fn values_along_x_works() {
        let mesh = Samples::one_tri6();
        let feat = Features::new(&mesh, false);
        let p1 = SampleParams::param_diffusion();
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut state = FemState::new(&input, &config).unwrap();
        state.uu[0] = 1.0;
        state.uu[1] = 2.0;
        state.uu[2] = 3.0;
        state.uu[3] = 4.0;
        state.uu[4] = 5.0;
        state.uu[5] = 6.0;
        let output = FemOutput::new(&input, None, None, None).unwrap();
        let (ids, xx, dd) = output.values_along_x(&feat, &state, Dof::T, 0.0, any_x).unwrap();
        assert_eq!(ids, &[0, 3, 1]);
        assert_eq!(xx, &[0.0, 0.5, 1.0]);
        assert_eq!(dd, &[1.0, 4.0, 2.0]);
    }
}
