use super::{FemBase, FileIo};
use crate::base::Dof;
use crate::fem::FemState;
use crate::StrError;
use gemlab::mesh::Mesh;
use std::fmt::Write;
use std::fs::File;
use std::io::Write as IoWrite;

impl FileIo {
    /// Writes a file associated with a single time station to perform visualization with ParaView
    ///
    /// The files will be indexed with `index` corresponding to each time station.
    pub fn write_vtu(&self, mesh: &Mesh, base: &FemBase, state: &FemState, index: usize) -> Result<(), StrError> {
        if !self.active {
            return Err("FileIo must be activated first");
        }

        let ndim = mesh.ndim;
        let npoint = mesh.points.len();
        let ncell = mesh.cells.len();
        if ncell < 1 {
            return Err("there are no cells to write");
        }

        // auxiliary information
        let enabled_dofs = base.dofs.enabled();
        let not_displacement_dof: Vec<_> = enabled_dofs
            .iter()
            .filter(|&&dof| !(dof == Dof::Ux || dof == Dof::Uy || dof == Dof::Uz))
            .copied()
            .collect();

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
        if enabled_dofs.contains(&Dof::Ux) {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n"
            )
            .unwrap();
            for point in &mesh.points {
                let ux = match base.dofs.eq(point.id, Dof::Ux).ok() {
                    Some(eq) => state.u[eq],
                    None => 0.0,
                };
                let uy = match base.dofs.eq(point.id, Dof::Uy).ok() {
                    Some(eq) => state.u[eq],
                    None => 0.0,
                };
                let uz = match base.dofs.eq(point.id, Dof::Uz).ok() {
                    Some(eq) => state.u[eq],
                    None => 0.0,
                };
                write!(&mut buffer, "{:?} {:?} {:?} ", ux, uy, uz).unwrap();
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }
        for dof in &not_displacement_dof {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"{:?}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                dof
            )
            .unwrap();
            for point in &mesh.points {
                let value = match base.dofs.eq(point.id, *dof).ok() {
                    Some(eq) => state.u[eq],
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

        // write file
        let path = self.path_vtu(index);
        let mut file = File::create(&path).map_err(|_| "cannot create VTU file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write VTU file")?;
        Ok(())
    }

    /// Writes a summary file for all time stations to perform visualization with ParaView
    pub fn write_pvd(&self) -> Result<(), StrError> {
        if !self.active {
            return Err("FileIo must be activated first");
        }

        // header
        let mut buffer = String::new();
        write!(&mut buffer, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n").unwrap();

        // add VTU entries to PVD file
        for index in &self.indices {
            let vtu_fn = self.path_vtu(*index);
            write!(
                &mut buffer,
                "<DataSet timestep=\"{:?}\" file=\"{}\" />\n",
                self.times[*index], vtu_fn
            )
            .unwrap();
        }

        // footer
        write!(&mut buffer, "</Collection>\n</VTKFile>\n").unwrap();

        // write file
        let path = self.path_pvd();
        let mut file = File::create(&path).map_err(|_| "cannot create PVD file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write PVD file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::{Config, Dof, Elem, Essential, ParamBeam, ParamPorousSldLiq, ParamSolid};
    use crate::fem::{FemBase, FemState, FileIo};
    use gemlab::mesh::Samples;
    use std::fs;

    #[test]
    fn write_vtu_captures_errors() {
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let file_io = FileIo::new();
        assert_eq!(
            file_io.write_vtu(&mesh, &base, &state, 0).err(),
            Some("FileIo must be activated first")
        );
    }

    #[test]
    fn write_vtu_works() {
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // Generates a displacement field corresponding to a simple shear deformation
        // Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
        let strain = 1.23;
        let npoint = mesh.points.len();
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            let eq = base.dofs.eq(p, Dof::Ux).unwrap();
            state.u[eq] = strain * y;
        }

        let index = 0;
        let fn_stem = "test_write_vtu_works";
        let mut file_io = FileIo::new();
        file_io.activate(&mesh, &base, "/tmp/pmsim", fn_stem).unwrap();
        file_io.write_vtu(&mesh, &base, &state, index).unwrap();

        let fn_path = file_io.path_vtu(index);
        let contents = fs::read_to_string(&fn_path).map_err(|_| "cannot open file").unwrap();
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

    #[test]
    fn write_vtu_works_mixed() {
        //                     {Ux→15}
        //    {Ux→21}          {Uy→16}
        //    {Uy→22}  {Ux→19} {Rz→17}
        //    {Pl→23}  {Uy→20} {Pl→18} {Ux→13}
        //         8------7------6._   {Uy→14}
        //         |       [3](3)|  '-.5
        //         |  [0]        |     '-._
        // {Ux→24} 9  (1)      *10  [1]    '4 {Ux→11}
        // {Uy→25} |             |  (2)  .-'  {Uy→12}
        //         |       [2](3)|   _.3'
        //         0------1------2.-'  {Ux→9}
        //     {Ux→0}  {Ux→3}  {Ux→5}  {Uy→10}
        //     {Uy→1}  {Uy→4}  {Uy→6}
        //     {Pl→2}          {Rz→7}
        //                     {Pl→8}
        //  *10 => {Ux→26, Uy→27, Rz→28}
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let base = FemBase::new(
            &mesh,
            [(1, Elem::PorousSldLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Beam(p3))],
        )
        .unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let mut state = FemState::new(&mesh, &base, &essential, &config).unwrap();

        // Generates a displacement field corresponding to a simple shear deformation
        // Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
        let strain = 1.23;
        let npoint = mesh.points.len();
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            let eq = base.dofs.eq(p, Dof::Ux).unwrap();
            state.u[eq] = strain * y;
        }

        // Applies liquid pressure proportional to the y coordinate
        for p in 0..npoint {
            let y = mesh.points[p].coords[1];
            if base.dofs.contains(p, Dof::Pl) {
                let eq = base.dofs.eq(p, Dof::Pl).unwrap();
                state.u[eq] = 100.0 * (1.0 + y);
            }
        }

        let index = 0;
        let fn_stem = "test_write_vtu_works_mixed";
        let mut file_io = FileIo::new();
        file_io.activate(&mesh, &base, "/tmp/pmsim", fn_stem).unwrap();
        file_io.write_vtu(&mesh, &base, &state, index).unwrap();

        let fn_path = file_io.path_vtu(index);
        let contents = fs::read_to_string(&fn_path).map_err(|_| "cannot open file").unwrap();
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">
<UnstructuredGrid>
<Piece NumberOfPoints="11" NumberOfCells="4">
<Points>
<DataArray type="Float64" NumberOfComponents="3" format="ascii">
0.0 0.0 0.0 0.5 0.0 0.0 1.0 0.0 0.0 1.433 0.25 0.0 1.866 0.5 0.0 1.433 0.75 0.0 1.0 1.0 0.0 0.5 1.0 0.0 0.0 1.0 0.0 0.0 0.5 0.0 1.0 0.5 0.0 
</DataArray>
</Points>
<Cells>
<DataArray type="Int32" Name="connectivity" format="ascii">
0 2 6 8 1 10 7 9 2 4 6 3 5 10 2 10 10 6 
</DataArray>
<DataArray type="Int32" Name="offsets" format="ascii">
8 14 16 18 
</DataArray>
<DataArray type="UInt8" Name="types" format="ascii">
23 22 3 3 
</DataArray>
</Cells>
<PointData Scalars="TheScalars">
<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii">
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.3075 0.0 0.0 0.615 0.0 0.0 0.9225 0.0 0.0 1.23 0.0 0.0 1.23 0.0 0.0 1.23 0.0 0.0 0.615 0.0 0.0 0.615 0.0 0.0 
</DataArray>
<DataArray type="Float64" Name="Rz" NumberOfComponents="1" format="ascii">
0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 
</DataArray>
<DataArray type="Float64" Name="Pl" NumberOfComponents="1" format="ascii">
100.0 0.0 100.0 0.0 0.0 0.0 200.0 0.0 200.0 0.0 0.0 
</DataArray>
</PointData>
</Piece>
</UnstructuredGrid>
</VTKFile>
"#
        );
    }

    #[test]
    fn write_pvd_captures_errors() {
        let file_io = FileIo::new();
        assert_eq!(file_io.write_pvd().err(), Some("FileIo must be activated first"));
    }

    #[test]
    fn write_pvd_works() {
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        let fn_stem = "test_write_pvd_works";
        let mut file_io = FileIo::new();

        file_io.activate(&mesh, &base, "/tmp/pmsim", fn_stem).unwrap();
        file_io.write_state(&state).unwrap();
        file_io.write_pvd().unwrap();

        let fn_path = file_io.path_pvd();
        let contents = fs::read_to_string(&fn_path).map_err(|_| "cannot open file").unwrap();
        assert_eq!(
            contents,
            r#"<?xml version="1.0"?>
<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">
<Collection>
<DataSet timestep="0.0" file="/tmp/pmsim/test_write_pvd_works-00000000000000000000.vtu" />
</Collection>
</VTKFile>
"#
        );
    }
}
