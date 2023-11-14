use crate::base::Dof;
use crate::fem::{Elements, FemInput, FemState};
use crate::StrError;
use gemlab::mesh::{At, Features, PointId};
use std::fs;

pub const DEFAULT_OUTPUT_DIR: &str = "/tmp/pmsim/results";

/// Assists in the post-processing of results
pub struct FemOutput<'a> {
    input: &'a FemInput<'a>,
    filename_stem: Option<String>,
    output_directory: String,
    output_count: usize,
    callback: Option<fn(&FemState, usize)>,
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
            None => DEFAULT_OUTPUT_DIR,
        };
        fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;

        // results
        Ok(FemOutput {
            input,
            filename_stem,
            output_directory: out_dir.to_string(),
            output_count: 0,
            callback,
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
            self.output_count += 1;
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
