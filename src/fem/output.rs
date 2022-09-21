use super::{Data, State};
use crate::base::Dof;
use crate::StrError;
use std::collections::HashSet;
use std::ffi::OsStr;
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

/// Output calculation results
pub struct Output<'a> {
    pub data: &'a Data<'a>,
    pub enabled_dofs: HashSet<Dof>,
    pub not_displacement_dof: Vec<Dof>,
}

impl<'a> Output<'a> {
    /// Allocate new instance
    pub fn new(data: &'a Data) -> Self {
        let mut enabled_dofs = HashSet::new();
        for map in &data.equations.all {
            for dof in map.keys() {
                enabled_dofs.insert(*dof);
            }
        }
        let not_displacement_dof: Vec<_> = enabled_dofs
            .iter()
            .filter(|&&dof| !(dof == Dof::Ux || dof == Dof::Uy || dof == Dof::Uz))
            .copied()
            .collect();
        Output {
            data,
            enabled_dofs,
            not_displacement_dof,
        }
    }

    /// Writes Paraview's VTU file
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_vtu<P>(&self, state: &State, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let mesh = self.data.mesh;
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
                let ux = match self.data.equations.eq(point.id, Dof::Ux).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                let uy = match self.data.equations.eq(point.id, Dof::Uy).ok() {
                    Some(eq) => state.uu[eq],
                    None => 0.0,
                };
                let uz = match self.data.equations.eq(point.id, Dof::Uz).ok() {
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
                let value = match self.data.equations.eq(point.id, *dof).ok() {
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

        // force sync
        file.sync_all().map_err(|_| "cannot sync file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Output;
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::{Data, State};
    use gemlab::mesh::Samples;
    use std::fs;

    #[test]
    fn new_works() {
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        let mut state = State::new(&data, &config).unwrap();

        // Generates a displacement field corresponding to a simple shear deformation
        // Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
        let strain = 1.23;
        let npoint = mesh.points.len();
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            state.uu[0 + mesh.ndim * p] = strain * y;
        }

        let output = Output::new(&data);
        let file_path = "/tmp/pmsim/test_output.vtu";
        output.write_vtu(&state, file_path).unwrap();

        let contents = fs::read_to_string(file_path).map_err(|_| "cannot open file").unwrap();
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
