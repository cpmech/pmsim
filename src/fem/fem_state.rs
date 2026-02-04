use super::SecondaryValues;
use crate::base::{Config, ElemType, Schema};
use crate::StrError;
use gemlab::integ::Gauss;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds the state of a simulation
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemState {
    /// Time t
    pub time: f64,

    /// Loading factor λ
    pub lambda: f64,

    /// Indicates whether a load reversal occurred
    pub reverse: bool,

    /// Holds the α1(time) coefficient for the dynamics method
    pub alpha1: f64,

    /// Holds the α2(time) coefficient for the dynamics method
    pub alpha2: f64,

    /// Holds the α3(time) coefficient for the dynamics method
    pub alpha3: f64,

    /// Holds the α4(time) coefficient for the dynamics method
    pub alpha4: f64,

    /// Holds the α5(time) coefficient for the dynamics method
    pub alpha5: f64,

    /// Holds the α6(time) coefficient for the dynamics method
    pub alpha6: f64,

    /// Holds the α7(time) coefficient for the dynamics method
    pub alpha7: f64,

    /// Holds the α8(time) coefficient for the dynamics method
    pub alpha8: f64,

    /// Holds the β1(time) coefficient for the transient/dynamics method
    pub beta1: f64,

    /// Holds the β2(time) coefficient for the dynamics method
    pub beta2: f64,

    /// Time increment Δt
    pub ddt: f64,

    /// Lambda increment Δλ
    pub ddl: f64,

    /// Primary variables increment ΔU
    ///
    /// (neq)
    pub dduu: Vector,

    /// Primary unknowns (using uu to indicate capital U)
    ///
    /// (neq)
    pub uu: Vector,

    /// First time derivative of primary unknowns dU/dt
    ///
    /// (neq)
    pub vv: Vector,

    /// Second time derivative of primary unknowns d²U/dt²
    ///
    /// (neq)
    pub aa: Vector,

    /// Auxiliary time-discretization variable (Theta method) U★
    ///
    /// (neq)
    pub uu_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method) V★
    ///
    /// (neq)
    pub vv_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method) A★
    ///
    /// (neq)
    pub aa_star: Vector,

    /// Holds the secondary values (e.g. stress) at all integration (Gauss) points of all elements
    ///
    /// (ncell)
    pub gauss: Vec<SecondaryValues>,
}

impl FemState {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, schema: &Schema, config: &Config) -> Result<FemState, StrError> {
        // check number of cells
        let ncell = mesh.cells.len();
        if ncell == 0 {
            return Err("there are no cells in the mesh");
        }

        // secondary values (e.g. stress) at all integration (Gauss) points of all elements
        let empty = SecondaryValues::new_empty();
        let mut gauss = vec![empty; ncell];

        // gather information about element types
        let mandel = config.ideal.mandel();
        let mut has_diffusion = false;
        let mut has_rod_or_beam = false;
        let mut has_solid = false;
        let mut has_porous_fluid = false;
        let mut has_porous_solid = false;
        for cell in &mesh.cells {
            let elem_type = schema.elem_type(cell.marker)?;
            let ngauss_opt = elem_type.ngauss();
            let ngauss = Gauss::new_or_sized(cell.kind, ngauss_opt)?.npoint();
            match elem_type {
                ElemType::Diffusion(..) => {
                    has_diffusion = true;
                    if config.model_settings(cell.marker).save_flux {
                        gauss[cell.id].allocate_diffusion(ngauss, mesh.ndim);
                    }
                }
                ElemType::Rod(..) => {
                    has_rod_or_beam = true;
                }
                ElemType::Beam(..) => {
                    has_rod_or_beam = true;
                }
                ElemType::Solid(param) => {
                    has_solid = true;
                    let n_int_var = param.n_int_var();
                    gauss[cell.id].allocate_solid(mandel, ngauss, n_int_var);
                }
                ElemType::PorousLiq(..) => {
                    has_porous_fluid = true;
                    gauss[cell.id].allocate_porous_liq(ngauss);
                }
                ElemType::PorousLiqGas(..) => {
                    has_porous_fluid = true;
                    gauss[cell.id].allocate_porous_liq_gas(ngauss);
                }
                ElemType::PorousSldLiq(param) => {
                    has_porous_solid = true;
                    let n_int_var = param.n_int_var();
                    gauss[cell.id].allocate_porous_sld_liq(mandel, ngauss, n_int_var);
                }
                ElemType::PorousSldLiqGas(param) => {
                    has_porous_solid = true;
                    let n_int_var = param.n_int_var();
                    gauss[cell.id].allocate_porous_sld_liq_gas(mandel, ngauss, n_int_var);
                }
            };
        }

        // check elements
        if has_diffusion && (has_rod_or_beam || has_solid || has_porous_fluid || has_porous_solid) {
            return Err("cannot combine Diffusion elements with other elements");
        }
        if has_porous_fluid && (has_diffusion || has_rod_or_beam || has_solid || has_porous_solid) {
            return Err("cannot combine PorousLiq or PorousLiqGas with other elements");
        }

        // number of equations = total number of DOFs
        let neq = schema.ndof()?;

        // primary variables
        let dduu = Vector::new(neq);
        let uu = Vector::new(neq);
        let (uu_star, vv, vv_star) = if config.transient || config.dynamics {
            (Vector::new(neq), Vector::new(neq), Vector::new(neq))
        } else {
            (Vector::new(0), Vector::new(0), Vector::new(0))
        };
        let (aa, aa_star) = if config.dynamics {
            (Vector::new(neq), Vector::new(neq))
        } else {
            (Vector::new(0), Vector::new(0))
        };

        // allocate new instance
        Ok(FemState {
            time: 0.0,
            lambda: 0.0,
            reverse: false,
            alpha1: 0.0,
            alpha2: 0.0,
            alpha3: 0.0,
            alpha4: 0.0,
            alpha5: 0.0,
            alpha6: 0.0,
            alpha7: 0.0,
            alpha8: 0.0,
            beta1: 0.0,
            beta2: 0.0,
            ddt: 0.0,
            ddl: 0.0,
            dduu,
            uu,
            vv,
            aa,
            uu_star,
            vv_star,
            aa_star,
            gauss,
        })
    }

    /// Reads a JSON file containing the state data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open state file")?;
        let buffered = BufReader::new(data);
        let state = serde_json::from_reader(buffered).map_err(|_| "cannot parse state file")?;
        Ok(state)
    }

    /// Writes a JSON file with the state data
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
        let mut file = File::create(&path).map_err(|_| "cannot create state file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write state file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemState;
    use crate::base::{new_empty_mesh_2d, Config, Schema};
    use crate::base::{ParamBeam, ParamDiffusion, ParamPorousLiq, ParamPorousLiqGas};
    use crate::base::{ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRod, ParamSolid};
    use gemlab::mesh::Samples;

    #[test]
    fn new_handles_errors() {
        let mesh = new_empty_mesh_2d();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &schema, &config).err(),
            Some("there are no cells in the mesh")
        );

        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamDiffusion::sample();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamRod::sample();
        let mut schema = Schema::new();
        schema
            .add_diffusion(1, p1)
            .add_solid(2, p2)
            .add_rod(3, p3)
            .build(&mesh)
            .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &schema, &config).err(),
            Some("cannot combine Diffusion elements with other elements")
        );

        let p1 = ParamPorousLiq::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema
            .add_porous_liq(1, p1)
            .add_solid(2, p2)
            .add_rod(3, p3)
            .build(&mesh)
            .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &schema, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );

        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema
            .add_porous_liq_gas(1, p1)
            .add_solid(2, p2)
            .add_rod(3, p3)
            .build(&mesh)
            .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &schema, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );
    }

    #[test]
    fn new_works_mixed() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let mut schema = Schema::new();
        schema
            .add_porous_sld_liq(1, p1)
            .add_solid(2, p2)
            .add_beam(3, p3)
            .build(&mesh)
            .unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
    }

    #[test]
    fn new_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let mut schema = Schema::new();
        schema.add_diffusion(1, p1).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.transient = true;
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
        assert_eq!(state.vv.dim(), schema.ndof().unwrap());
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), schema.ndof().unwrap());
        assert_eq!(state.vv_star.dim(), schema.ndof().unwrap());
        assert_eq!(state.aa_star.dim(), 0);
    }

    #[test]
    fn new_works_rod_only() {
        let mesh = Samples::one_lin2();
        let p1 = ParamRod::sample();
        let mut schema = Schema::new();
        schema.add_rod(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
        assert_eq!(state.vv.dim(), 0);
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), 0);
        assert_eq!(state.vv_star.dim(), 0);
        assert_eq!(state.aa_star.dim(), 0);
    }

    #[test]
    fn new_works_porous_liq() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousLiq::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema.add_porous_liq(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
    }

    #[test]
    fn new_works_porous_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let mut schema = Schema::new();
        schema.add_porous_liq_gas(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
    }

    #[test]
    fn new_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let mut schema = Schema::new();
        schema.add_porous_sld_liq_gas(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = ParamRod::sample();
        let p2 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_rod(1, p1).add_solid(2, p2).build(&mesh).unwrap();
        let mut config = Config::new(&mesh);
        config.dynamics = true;
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        assert_eq!(state.dduu.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu.dim(), schema.ndof().unwrap());
        assert_eq!(state.vv.dim(), schema.ndof().unwrap());
        assert_eq!(state.aa.dim(), schema.ndof().unwrap());
        assert_eq!(state.uu_star.dim(), schema.ndof().unwrap());
        assert_eq!(state.vv_star.dim(), schema.ndof().unwrap());
        assert_eq!(state.aa_star.dim(), schema.ndof().unwrap());
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_lin2();
        let p1 = ParamRod::sample();
        let mut schema = Schema::new();
        schema.add_rod(1, p1).build(&mesh).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &schema, &config).unwrap();
        let clone = state.clone();
        let str_ori = format!("{:?}", clone).to_string();
        assert_eq!(format!("{:?}", clone), str_ori);
        // serialize
        let json = serde_json::to_string(&clone).unwrap();
        // deserialize
        let read: FemState = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), str_ori);
    }
}
