use super::{FemInput, SecondaryValues};
use crate::base::{Config, Element};
use crate::material::LocalState;
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::Vector;
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds state of a simulation, including primary and secondary variables
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

    /// Secondary values at all elements and all integration points (output only)
    ///
    /// (n_cells)
    pub secondary_values: Option<Vec<SecondaryValues>>,
}

impl FemState {
    /// Allocates a new instance
    pub fn new(input: &FemInput, config: &Config) -> Result<FemState, StrError> {
        // check number of cells
        if input.mesh.cells.len() == 0 {
            return Err("there are no cells in the mesh");
        }

        // gather information about element types
        let mut has_diffusion = false;
        let mut has_rod_or_beam = false;
        let mut has_solid = false;
        let mut has_porous_fluid = false;
        let mut has_porous_solid = false;
        for cell in &input.mesh.cells {
            let element = input.attributes.get(cell).unwrap(); // already checked by Data
            match element {
                Element::Diffusion(..) => has_diffusion = true,
                Element::Rod(..) => has_rod_or_beam = true,
                Element::Beam(..) => has_rod_or_beam = true,
                Element::Solid(..) => has_solid = true,
                Element::PorousLiq(..) => has_porous_fluid = true,
                Element::PorousLiqGas(..) => has_porous_fluid = true,
                Element::PorousSldLiq(..) => has_porous_solid = true,
                Element::PorousSldLiqGas(..) => has_porous_solid = true,
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
        let n_equation = input.equations.n_equation;

        // primary variables
        let t = config.control.t_ini;
        let dt = (config.control.dt)(t);
        let duu = Vector::new(n_equation);
        let uu = Vector::new(n_equation);
        let (uu_star, vv, vv_star) = if config.transient || config.dynamics {
            (
                Vector::new(n_equation),
                Vector::new(n_equation),
                Vector::new(n_equation),
            )
        } else {
            (Vector::new(0), Vector::new(0), Vector::new(0))
        };
        let (aa, aa_star) = if config.dynamics {
            (Vector::new(n_equation), Vector::new(n_equation))
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
            secondary_values: None,
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
        let input = File::open(path).map_err(|_| "cannot open file")?;
        let buffered = BufReader::new(input);
        let state = serde_json::from_reader(buffered).map_err(|_| "cannot parse JSON file")?;
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
        let mut file = File::create(&path).map_err(|_| "cannot create file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write file")?;
        Ok(())
    }

    /// Extracts stresses, strains, and internal values
    ///
    /// Returns `(stress_state, epsilon)`
    ///
    /// The results are only available if [crate::fem::FemOutput] was provided a `filename_stem`
    /// and [Config::out_secondary_values] was set true.
    ///
    /// **Note:** This function will return an error if stresses/strains are not available (e.g. in a Diffusion simulation)
    /// or not present in the integration (Gauss) point of the selected cell/element.
    pub fn extract_stresses_and_strains(
        &self,
        cell_id: CellId,
        integ_point: usize,
    ) -> Result<&LocalState, StrError> {
        match &self.secondary_values {
            Some(all_values) => {
                let values = &all_values[cell_id];
                let state = match &values.stresses_and_strains {
                    Some(states) => &states.all[integ_point],
                    None => return Err("element does not have stresses at the selected integration point"),
                };
                Ok(state)
            }
            None => {
                Err("secondary values are not available for this problem (need config.out_secondary_values = true)")
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemState;
    use crate::base::{new_empty_mesh_2d, Config, Element, SampleParams};
    use crate::fem::FemInput;
    use gemlab::mesh::Samples;

    #[test]
    fn new_handles_errors() {
        let empty_mesh = new_empty_mesh_2d();
        let p1 = SampleParams::param_solid();
        let input = FemInput::new(&empty_mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new(&empty_mesh);
        assert_eq!(
            FemState::new(&input, &config).err(),
            Some("there are no cells in the mesh")
        );

        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_diffusion();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_rod();
        let input = FemInput::new(
            &mesh,
            [
                (1, Element::Diffusion(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&input, &config).err(),
            Some("cannot combine Diffusion elements with other elements")
        );

        let p1 = SampleParams::param_porous_liq();
        let input = FemInput::new(
            &mesh,
            [
                (1, Element::PorousLiq(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&input, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );

        let p1 = SampleParams::param_porous_liq_gas();
        let input = FemInput::new(
            &mesh,
            [
                (1, Element::PorousLiqGas(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new(&mesh);
        assert_eq!(
            FemState::new(&input, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );
    }

    #[test]
    fn new_works_mixed() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid_drucker_prager();
        let p3 = SampleParams::param_beam();
        let input = FemInput::new(
            &mesh,
            [
                (1, Element::PorousSldLiq(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Beam(p3)),
            ],
        )
        .unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.dt, 1.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
    }

    #[test]
    fn new_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let mut config = Config::new(&mesh);
        config.transient = true;
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
        assert_eq!(state.vv.dim(), input.equations.n_equation);
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), input.equations.n_equation);
        assert_eq!(state.vv_star.dim(), input.equations.n_equation);
        assert_eq!(state.aa_star.dim(), 0);
    }

    #[test]
    fn new_works_rod_only() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let input = FemInput::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
        assert_eq!(state.vv.dim(), 0);
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), 0);
        assert_eq!(state.vv_star.dim(), 0);
        assert_eq!(state.aa_star.dim(), 0);
    }

    #[test]
    fn new_works_porous_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq();
        let input = FemInput::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
    }

    #[test]
    fn new_works_porous_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq_gas();
        let input = FemInput::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
    }

    #[test]
    fn new_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let input = FemInput::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = SampleParams::param_rod();
        let p2 = SampleParams::param_solid();
        let input = FemInput::new(&mesh, [(1, Element::Rod(p1)), (2, Element::Solid(p2))]).unwrap();
        let mut config = Config::new(&mesh);
        config.dynamics = true;
        let state = FemState::new(&input, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.duu.dim(), input.equations.n_equation);
        assert_eq!(state.uu.dim(), input.equations.n_equation);
        assert_eq!(state.vv.dim(), input.equations.n_equation);
        assert_eq!(state.aa.dim(), input.equations.n_equation);
        assert_eq!(state.uu_star.dim(), input.equations.n_equation);
        assert_eq!(state.vv_star.dim(), input.equations.n_equation);
        assert_eq!(state.aa_star.dim(), input.equations.n_equation);
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let input = FemInput::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new(&mesh);
        let state_ori = FemState::new(&input, &config).unwrap();
        let state = state_ori.clone();
        let str_ori = format!("{:?}", state).to_string();
        assert!(str_ori.len() > 0);
        // serialize
        let json = serde_json::to_string(&state).unwrap();
        // deserialize
        let read: FemState = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), str_ori);
    }
}
