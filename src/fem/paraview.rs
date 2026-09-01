use crate::base::Dof;
use crate::fem::{FemState, PostProc, PostProcMemo};
use crate::StrError;
use std::fmt::Write;
use std::fs::{self, File};
use std::io::Write as IoWrite;
use std::path::Path;

const WITH_CELL_DATA: bool = false;

impl<const DIM: usize> PostProc<DIM> {
    /// Writes a file associated with a single time station to perform visualization with ParaView
    ///
    /// **Warning:** This function **does not** create the output directory if it does not exist.
    ///
    /// Returns the path to the VTU file
    ///
    /// The files will be indexed with `index` corresponding to each time station.
    pub fn write_vtu(
        &self,
        memo: &mut PostProcMemo,
        dir: &str,
        fn_stem: &str,
        state: &FemState<DIM>,
        index: usize,
        with_elastic_flags: bool,
    ) -> Result<String, StrError> {
        // auxiliary variables
        let ndim = self.mesh.ndim;
        let npoint = self.mesh.points.len();
        let ncell = self.mesh.cells.len();
        if ncell < 1 {
            return Err("there are no cells to write");
        }

        // DOF information
        let (displacement_dofs, non_displacement_dofs) = self.schema.enabled_dofs()?;

        // detect features
        let mut has_phi_flux = false;
        let mut has_elastic_flag = false;
        for g in &state.gauss {
            if g.diffusion.len() > 0 {
                has_phi_flux = true;
            }
            if g.solid.len() > 0 || g.porous_sld_liq.len() > 0 || g.porous_sld_liq_gas.len() > 0 {
                has_elastic_flag = true;
            }
        }

        // extrapolate flux from Gauss points to points
        let ww_at_nodes = if has_phi_flux {
            let cell_ids: Vec<_> = (0..ncell).into_iter().collect();
            Some(self.nodal_fluxes_patch(memo, state, &cell_ids, Dof::Phi, |_, _, _| true)?)
        } else {
            None
        };

        // elastic flags
        let elastic_flags = if with_elastic_flags && has_elastic_flag {
            let cell_ids: Vec<_> = (0..ncell).into_iter().collect();
            Some(self.gauss_elastic_flags_patch(memo, state, &cell_ids, |_, _, _| true)?)
        } else {
            None
        };

        // output buffer
        let mut buffer = String::new();

        // header ----------------------------------------------------------------------------------------------------------------
        let (npoint_total, ncell_total) = if elastic_flags.is_some() {
            let ngauss = elastic_flags.as_ref().unwrap().xx.len(); // each Gauss point corresponds to a point and a single-vertex cell
            (npoint + ngauss, ncell + ngauss)
        } else {
            (npoint, ncell)
        };
        write!(
            &mut buffer,
            "<?xml version=\"1.0\"?>\n\
             <VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n\
             <UnstructuredGrid>\n\
             <Piece NumberOfPoints=\"{}\" NumberOfCells=\"{}\">\n",
            npoint_total, ncell_total
        )
        .unwrap();

        // topology --------------------------------------------------------------------------------------------------------------

        // nodes: coordinates
        write!(
            &mut buffer,
            "<Points>\n\
             <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">\n",
        )
        .unwrap();
        for index in 0..npoint {
            for dim in 0..ndim {
                write!(&mut buffer, "{:?} ", self.mesh.points[index].coords[dim]).unwrap();
            }
            if ndim == 2 {
                write!(&mut buffer, "0.0 ").unwrap();
            }
        }
        if let Some(flags) = elastic_flags.as_ref() {
            for k in 0..flags.k_to_id.len() {
                write!(&mut buffer, "{:?} {:?} ", flags.xx[k], flags.yy[k]).unwrap();
                if ndim == 2 {
                    write!(&mut buffer, "0.0 ").unwrap();
                } else {
                    write!(&mut buffer, "{:?} ", flags.zz[k]).unwrap();
                }
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
        for cell in &self.mesh.cells {
            if cell.kind.vtk_type().is_none() {
                return Err("cannot generate VTU file because VTK cell type is not available");
            }
            for p in &cell.points {
                write!(&mut buffer, "{} ", p).unwrap();
            }
        }
        if let Some(flags) = elastic_flags.as_ref() {
            for k in 0..flags.k_to_id.len() {
                write!(&mut buffer, "{} ", npoint + k).unwrap();
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
        for cell in &self.mesh.cells {
            offset += cell.points.len();
            write!(&mut buffer, "{} ", offset).unwrap();
        }
        if let Some(flags) = elastic_flags.as_ref() {
            for _ in 0..flags.k_to_id.len() {
                offset += 1; // single-vertex cell
                write!(&mut buffer, "{} ", offset).unwrap();
            }
        }

        // elements: types
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n"
        )
        .unwrap();
        for cell in &self.mesh.cells {
            if let Some(vtk) = cell.kind.vtk_type() {
                write!(&mut buffer, "{} ", vtk).unwrap();
            }
        }
        if let Some(flags) = elastic_flags.as_ref() {
            for _ in 0..flags.k_to_id.len() {
                write!(&mut buffer, "1 ").unwrap(); // VTK_VERTEX = 1
            }
        }
        write!(
            &mut buffer,
            "\n</DataArray>\n\
             </Cells>\n"
        )
        .unwrap();

        // begin point data section ----------------------------------------------------------------------------------------------
        write!(&mut buffer, "<PointData Scalars=\"TheScalars\">\n").unwrap();

        // point data: displacement DOFs
        if !displacement_dofs.is_empty() {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"displacement\" NumberOfComponents=\"3\" format=\"ascii\">\n"
            )
            .unwrap();
            for point in &self.mesh.points {
                let ux = if self.schema.has_dof(point.id, Dof::Ux)? {
                    let i = self.schema.dof_number(point.id, Dof::Ux)?;
                    state.uu[i]
                } else {
                    0.0
                };
                let uy = if self.schema.has_dof(point.id, Dof::Uy)? {
                    let i = self.schema.dof_number(point.id, Dof::Uy)?;
                    state.uu[i]
                } else {
                    0.0
                };
                let uz = if self.schema.has_dof(point.id, Dof::Uz)? {
                    let i = self.schema.dof_number(point.id, Dof::Uz)?;
                    state.uu[i]
                } else {
                    0.0
                };
                write!(&mut buffer, "{:?} {:?} {:?} ", ux, uy, uz).unwrap();
            }
            if let Some(flags) = elastic_flags.as_ref() {
                for _ in 0..flags.k_to_id.len() {
                    write!(&mut buffer, "0.0 0.0 0.0 ").unwrap(); // UNAVAILABLE
                }
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }

        // point data: non-displacement DOFs
        for dof in non_displacement_dofs {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"{:?}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                dof
            )
            .unwrap();
            for point in &self.mesh.points {
                let value = if self.schema.has_dof(point.id, *dof)? {
                    let i = self.schema.dof_number(point.id, *dof)?;
                    state.uu[i]
                } else {
                    0.0
                };
                write!(&mut buffer, "{:?} ", value).unwrap();
            }
            if let Some(flags) = elastic_flags.as_ref() {
                for _ in 0..flags.k_to_id.len() {
                    write!(&mut buffer, "0.0 0.0 0.0 ").unwrap(); // UNAVAILABLE
                }
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }

        // point data: flow vectors @ nodes
        if let Some(data) = ww_at_nodes.as_ref() {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"3\" format=\"ascii\">\n",
                data.label
            )
            .unwrap();
            for point in &self.mesh.points {
                let k = data.id_to_k[&point.id];
                let vx = data.vvx[k];
                let vy = data.vvy[k];
                let vz = if ndim == 3 { data.vvz[k] } else { 0.0 };
                write!(&mut buffer, "{:?} {:?} {:?} ", vx, vy, vz).unwrap();
            }
            if let Some(flags) = elastic_flags.as_ref() {
                for _ in 0..flags.k_to_id.len() {
                    write!(&mut buffer, "0.0 0.0 0.0 ").unwrap(); // UNAVAILABLE
                }
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }

        // point data: elastic flags @ the new vertices created from the Gauss points
        if let Some(flags) = elastic_flags.as_ref() {
            write!(
                &mut buffer,
                "<DataArray type=\"Float64\" Name=\"{}\" NumberOfComponents=\"1\" format=\"ascii\">\n",
                flags.label
            )
            .unwrap();
            for _ in &self.mesh.points {
                write!(&mut buffer, "-1.0 ").unwrap(); // UNAVAILABLE
            }
            for k in 0..flags.k_to_id.len() {
                write!(&mut buffer, "{:?} ", flags.values[k]).unwrap();
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
        }

        // end point data section ------------------------------------------------------------------------------------------------
        write!(&mut buffer, "</PointData>\n").unwrap();

        // cell data -------------------------------------------------------------------------------------------------------------
        if WITH_CELL_DATA {
            write!(&mut buffer, "<CellData Scalars=\"TheScalars\">\n").unwrap();

            // cell data: cell IDs
            write!(
                &mut buffer,
                "<DataArray type=\"Int32\" Name=\"cell_id\" NumberOfComponents=\"1\" format=\"ascii\">\n"
            )
            .unwrap();
            for cell in &self.mesh.cells {
                write!(&mut buffer, "{:?} ", cell.id).unwrap();
            }
            if let Some(flags) = elastic_flags.as_ref() {
                for k in 0..flags.k_to_id.len() {
                    write!(&mut buffer, "{} ", flags.k_to_id[k]).unwrap();
                }
            }
            write!(&mut buffer, "\n</DataArray>\n").unwrap();
            write!(&mut buffer, "</CellData>\n").unwrap();
        }

        // footer ----------------------------------------------------------------------------------------------------------------
        write!(
            &mut buffer,
            "</Piece>\n\
             </UnstructuredGrid>\n\
             </VTKFile>\n"
        )
        .unwrap();

        // write file
        let path = format!("{}/{}-{}.vtu", dir, fn_stem, index);
        let mut file = File::create(&path).map_err(|_| "cannot create VTU file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write VTU file")?;
        Ok(path)
    }

    /// Writes a summary file for all time stations to perform visualization with ParaView
    ///
    /// **Note:** This function creates the output directory if it does not exist.
    ///
    /// Returns the path to the PVD file
    // pub fn write_pvd(dir: &str, fn_stem: &str, indices: &[usize], times: &[f64]) -> Result<String, StrError> {
    pub fn write_pvd(&self, dir: &str, fn_stem: &str) -> Result<String, StrError> {
        // get indices and times
        let indices: Vec<_> = (0..self.files.nfile()).into_iter().collect();
        let times = &self.files.stations();

        // header
        let mut buffer = String::new();
        write!(&mut buffer, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n").unwrap();

        // add VTU entries to PVD file
        for index in indices {
            let vtu_fn = format!("{}/{}-{}.vtu", dir, fn_stem, index);
            write!(
                &mut buffer,
                "<DataSet timestep=\"{:?}\" file=\"{}\" />\n",
                times[index], vtu_fn
            )
            .unwrap();
        }

        // footer
        write!(&mut buffer, "</Collection>\n</VTKFile>\n").unwrap();

        // create directory and write file
        let full_path = format!("{}/{}.pvd", dir, fn_stem);
        let path = Path::new(&full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create PVD file")?;
        file.write_all(buffer.as_bytes()).map_err(|_| "cannot write PVD file")?;
        Ok(full_path)
    }

    /// Loads all states and writes Paraview's VTU and PVD files
    ///
    /// Returns the path to the PVD file
    pub fn write_paraview(
        &self,
        memo: &mut PostProcMemo,
        dir: &str,
        fn_stem: &str,
        with_elastic_flags: bool,
    ) -> Result<String, StrError> {
        // write VTU files
        for index in 0..self.nfile() {
            let state = self.read_file(index)?;
            self.write_vtu(memo, dir, fn_stem, &state, index, with_elastic_flags)?;
        }

        // write PVD file
        self.write_pvd(dir, fn_stem)
    }
}
