use crate::base::Dof;
use crate::fem::{Elements, FemInput, FemState};
use crate::StrError;
use gemlab::mesh::{At, Features, PointId};
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

pub const DEFAULT_OUTPUT_DIR: &str = "/tmp/pmsim/results";

/// Assists in the post-processing of results
pub struct FemOutput<'a> {
    input: &'a FemInput<'a>,
    enabled_dofs: HashSet<Dof>,
    not_displacement_dof: Vec<Dof>,
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
        // auxiliary information
        let mut enabled_dofs = HashSet::new();
        for map in &input.equations.all {
            for dof in map.keys() {
                enabled_dofs.insert(*dof);
            }
        }
        let not_displacement_dof: Vec<_> = enabled_dofs
            .iter()
            .filter(|&&dof| !(dof == Dof::Ux || dof == Dof::Uy || dof == Dof::Uz))
            .copied()
            .collect();

        // output directory
        let out_dir = match output_directory {
            Some(d) => d,
            None => DEFAULT_OUTPUT_DIR,
        };
        fs::create_dir_all(out_dir).map_err(|_| "cannot create output directory")?;

        // results
        Ok(FemOutput {
            input,
            enabled_dofs,
            not_displacement_dof,
            filename_stem,
            output_directory: out_dir.to_string(),
            output_count: 0,
            callback,
        })
    }

    /// Writes the current FEM state to a binary file and the results to a VTU file
    ///
    /// **Note:** No output is generated if `filename_stem` is None.
    pub(crate) fn write(&mut self, state: &mut FemState, elements: &mut Elements) -> Result<(), StrError> {
        if let Some(fn_stem) = &self.filename_stem {
            elements.output_internal_values(state)?;
            let path_state = format!("{}/{}-{:0>20}.state", self.output_directory, fn_stem, self.output_count);
            let path_vtu = format!("{}/{}-{:0>20}.vtu", self.output_directory, fn_stem, self.output_count);
            state.write(&path_state)?;
            self.write_vtu(state, &path_vtu)?;
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

    /// Writes Paraview's VTU file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub(crate) fn write_vtu<P>(&self, state: &FemState, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = self.input.mesh;
        let ndim = mesh.ndim;
        let npoint = mesh.points.len();
        let ncell = mesh.cells.len();
        if ncell < 1 {
            return Err("there are no cells to write");
        }

        // output buffer
        let mut buffer = String::new();

        // header
        write!(
            &mut buffer,
            "<?xml version=\"1.0\"?>\n\
             <VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
             <UnstructuredGrid>\n\
             <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            npoint, ncell
        )
        .unwrap();

        // nodes: coordinates
        write!(
            &mut buffer,
            "<Points>\n\
             <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        )
        .unwrap();
        for index in 0..npoint {
            for dim in 0..ndim {
                write!(&mut buffer, "{:?} ", mesh.points[index].coords[dim]).unwrap();
            }
            if ndim == 2 {
                write!(&mut buffer, "0.0 ").unwrap();
            }
        }
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             </Points>\n"
        )
        .unwrap();

        // elements: connectivity
        write!(
            &mut buffer,
            "<Cells>\n\
             <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">\n"
        )
        .unwrap();
        for cell in &mesh.cells {
            if cell.kind.vtk_type().is_none() {
                return Err("cannot generate VTU file because VTK cell type is not available");
            }
            for p in &cell.points {
                write!(&mut buffer, "{} ", p).unwrap();
            }
        }

        // elements: offsets
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">\n"
        )
        .unwrap();
        let mut offset = 0;
        for cell in &mesh.cells {
            offset += cell.points.len();
            write!(&mut buffer, "{} ", offset).unwrap();
        }

        // elements: types
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
        )
        .unwrap();
        for cell in &mesh.cells {
            if let Some(vtk) = cell.kind.vtk_type() {
                write!(&mut buffer, "{} ", vtk).unwrap();
            }
        }
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             </Cells>\n"
        )
        .unwrap();

        // data: points
        write!(&mut buffer, "<PointData Scalars=\"TheScalars\">\n").unwrap();
        if self.enabled_dofs.contains(&Dof::Ux) {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n"
            )
            .unwrap();
            for point in &mesh.points {
                let ux = match self.input.equations.eq(point.id, Dof::Ux).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                let uy = match self.input.equations.eq(point.id, Dof::Uy).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                let uz = match self.input.equations.eq(point.id, Dof::Uz).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                write!(&mut buffer, "{:?} {:?} {:?} ", ux, uy, uz).unwrap();
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }
        for dof in &self.not_displacement_dof {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"{:?}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                dof
            )
            .unwrap();
            for point in &mesh.points {
                let value = match self.input.equations.eq(point.id, *dof).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                write!(&mut buffer, "{:?} ", value).unwrap();
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }
        write!(&mut buffer, "</PointData>\n").unwrap();

        // footer
        write!(
            &mut buffer,
            "</Piece>\n\
             </UnstructuredGrid>\n\
             </VTKFile>\n"
        )
        .unwrap();

        // create directory
        let path = Path::new(full_path);
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }

        // write file
        let mut file = File::create(path).map_err(|_| "cannot create file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write file")?;
        file.sync_all().map_err(|_| "cannot sync file")
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
    use std::fs;

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

    #[test]
    fn write_vtu_works() {
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let mut state = FemState::new(&input, &config).unwrap();

        // Generates a displacement field corresponding to a simple shear deformation
        // Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
        let strain = 1.23;
        let npoint = mesh.points.len();
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            state.uu[0 + mesh.ndim * p] = strain * y;
        }

        let output = FemOutput::new(&input, Some("test_write_vtu_works".to_string()), None, None).unwrap();
        let path = format!(
            "{}/{}-{:0>20}.vtu",
            output.output_directory,
            output.filename_stem.as_ref().unwrap(),
            output.output_count
        );
        output.write_vtu(&state, &path).unwrap();

        let contents = fs::read_to_string(path).map_err(|_| "cannot open file").unwrap();
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="5" NumberOfCells="3">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0.0 0.2 0.0 1.2 0.0 0.0 2.2 0.1 0.0 1.8 1.0 0.0 0.5 1.2 0.0 
</DataArray>
</Points>
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
0 1 4 1 3 4 1 2 3 
</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">
3 6 9 
</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">
5 5 5 
</DataArray>
</Cells>
<PointData Scalars="TheScalars">
<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii">
0.246 0.0 0.0 0.0 0.0 0.0 0.123 0.0 0.0 1.23 0.0 0.0 1.476 0.0 0.0 
</DataArray>
</PointData>
</Piece>
</UnstructuredGrid>
</VTKFile>
"#
        );
    }
}
