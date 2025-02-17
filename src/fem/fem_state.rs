use super::{FemBase, SecondaryValues};
use crate::base::{Config, Elem, Essential};
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
    /// Time
    pub t: f64,

    /// Delta time
    pub dt: f64,

    /// Cumulated (for one timestep) primary unknowns {ΔU}
    pub duu: Vector,

    /// Primary unknowns {U}
    ///
    /// (n_equation)
    pub uu: Vector,

    /// First time derivative of primary unknowns d{U}/dt
    ///
    /// (n_equation)
    pub vv: Vector,

    /// Second time derivative of primary unknowns d²{U}/dt²
    ///
    /// (n_equation)
    pub aa: Vector,

    /// Auxiliary time-discretization variable (Theta method)
    ///
    /// (n_equation)
    pub uu_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method)
    ///
    /// (n_equation)
    pub vv_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method)
    ///
    /// (n_equation)
    pub aa_star: Vector,

    /// Holds the secondary values (e.g. stress) at all integration (Gauss) points of all elements
    ///
    /// (ncell)
    pub gauss: Vec<SecondaryValues>,
}

impl FemState {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, base: &FemBase, essential: &Essential, config: &Config) -> Result<FemState, StrError> {
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
            let elem = base.amap.get(cell.attribute).unwrap(); // already checked by Data
            let ngauss_opt = base.amap.ngauss(cell.attribute).unwrap();
            let ngauss = Gauss::new_or_sized(cell.kind, ngauss_opt)?.npoint();
            match elem {
                Elem::Diffusion(..) => {
                    has_diffusion = true;
                }
                Elem::Rod(..) => {
                    has_rod_or_beam = true;
                }
                Elem::Beam(..) => {
                    has_rod_or_beam = true;
                }
                Elem::Solid(param) => {
                    has_solid = true;
                    let n_int_var = param.n_int_var();
                    gauss[cell.id].allocate_solid(mandel, ngauss, n_int_var);
                }
                Elem::PorousLiq(..) => {
                    has_porous_fluid = true;
                    gauss[cell.id].allocate_porous_liq(ngauss);
                }
                Elem::PorousLiqGas(..) => {
                    has_porous_fluid = true;
                    gauss[cell.id].allocate_porous_liq_gas(ngauss);
                }
                Elem::PorousSldLiq(param) => {
                    has_porous_solid = true;
                    let n_int_var = param.n_int_var();
                    gauss[cell.id].allocate_porous_sld_liq(mandel, ngauss, n_int_var);
                }
                Elem::PorousSldLiqGas(param) => {
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

        // constants
        let neq_total = if config.lagrange_mult_method {
            let n_lagrange = essential.n_prescribed();
            base.dofs.size() + n_lagrange
        } else {
            base.dofs.size()
        };

        // primary variables
        let t = config.t_ini;
        let dt = (config.dt)(t);
        let duu = Vector::new(neq_total);
        let uu = Vector::new(neq_total);
        let (uu_star, vv, vv_star) = if config.transient || config.dynamics {
            (Vector::new(neq_total), Vector::new(neq_total), Vector::new(neq_total))
        } else {
            (Vector::new(0), Vector::new(0), Vector::new(0))
        };
        let (aa, aa_star) = if config.dynamics {
            (Vector::new(neq_total), Vector::new(neq_total))
        } else {
            (Vector::new(0), Vector::new(0))
        };

        // allocate new instance
        return Ok(FemState {
            t,
            dt,
            duu,
            uu,
            vv,
            aa,
            uu_star,
            vv_star,
            aa_star,
            gauss,
        });
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
    use crate::base::{new_empty_mesh_2d, Config, Elem, Essential};
    use crate::base::{ParamBeam, ParamDiffusion, ParamPorousLiq, ParamPorousLiqGas};
    use crate::base::{ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRod, ParamSolid};
    use crate::fem::FemBase;
    use gemlab::mesh::Samples;

    #[test]
    fn new_handles_errors() {
        let mesh = new_empty_mesh_2d();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &base, &essential, &config).err(),
            Some("there are no cells in the mesh")
        );

        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamDiffusion::sample();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamRod::sample();
        let base = FemBase::new(
            &mesh,
            [(1, Elem::Diffusion(p1)), (2, Elem::Solid(p2)), (3, Elem::Rod(p3))],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &base, &essential, &config).err(),
            Some("cannot combine Diffusion elements with other elements")
        );

        let p1 = ParamPorousLiq::sample_brooks_corey_constant();
        let base = FemBase::new(
            &mesh,
            [(1, Elem::PorousLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Rod(p3))],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &base, &essential, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );

        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let base = FemBase::new(
            &mesh,
            [(1, Elem::PorousLiqGas(p1)), (2, Elem::Solid(p2)), (3, Elem::Rod(p3))],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&mesh, &base, &essential, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );
    }

    #[test]
    fn new_works_mixed() {
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
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.dt, 1.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
    }

    #[test]
    fn new_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let essential = Essential::new();
        let mut config = Config::new(&mesh);
        config.transient = true;
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
        assert_eq!(state.vv.dim(), base.dofs.size());
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), base.dofs.size());
        assert_eq!(state.vv_star.dim(), base.dofs.size());
        assert_eq!(state.aa_star.dim(), 0);
    }

    #[test]
    fn new_works_rod_only() {
        let mesh = Samples::one_lin2();
        let p1 = ParamRod::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
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
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiq(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
    }

    #[test]
    fn new_works_porous_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousLiqGas::sample_brooks_corey_constant();
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiqGas(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
    }

    #[test]
    fn new_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::PorousSldLiqGas(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = ParamRod::sample();
        let p2 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Rod(p1)), (2, Elem::Solid(p2))]).unwrap();
        let essential = Essential::new();
        let mut config = Config::new(&mesh);
        config.dynamics = true;
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), base.dofs.size());
        assert_eq!(state.uu.dim(), base.dofs.size());
        assert_eq!(state.vv.dim(), base.dofs.size());
        assert_eq!(state.aa.dim(), base.dofs.size());
        assert_eq!(state.uu_star.dim(), base.dofs.size());
        assert_eq!(state.vv_star.dim(), base.dofs.size());
        assert_eq!(state.aa_star.dim(), base.dofs.size());
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_lin2();
        let p1 = ParamRod::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Rod(p1))]).unwrap();
        let essential = Essential::new();
        let config = Config::new(&mesh);
        let state = FemState::new(&mesh, &base, &essential, &config).unwrap();
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
