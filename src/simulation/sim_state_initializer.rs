#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use super::{Configuration, StateStress};
use crate::StrError;
use russell_lab::Vector;

pub enum IniOption {
    /// Geostatic initial state
    ///
    /// # Note
    ///
    /// * The valued parameter is an overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial isotropic stress state with σ_xx = σ_yy = σ_zz = value
    IsotropicStress,

    /// Zero initial state
    Zero,
}

pub struct SimStateInitializer {
    // config: &'a SimConfig<'a>, // Access to configuration
// geostatics: Option<Geostatics>,
// model_water: Option<ModelRealDensity>,
}

impl SimStateInitializer {
    pub fn new(config: &Configuration) -> Result<Self, StrError> {
        // let model_water = match &config.ini_option {
        //     IniOption::Geostatic(param) => Some(ModelRealDensity::new(param)?),
        //     _ => None,
        // };
        // Ok(SimStateInitializer { config, model_water })
        Err("TODO")
    }

    pub fn initialize_stress(&self, state: &mut StateStress, ip_coords: &Vector) -> Result<(), StrError> {
        /*
        match self.config.ini_option {
            IniOption::Geostatic(..) => {
                if let Some(model_water) = &self.model_water {
                    let space_ndim = self.config.mesh.space_ndim;
                    let elevation = ip_coords[space_ndim - 1];
                    let height = self.config.mesh.max[space_ndim - 1];
                    let gravity = self.config.gravity;
                    let pl = model_water.pressure_at_elevation(elevation, height, gravity)?;
                }
            }
            IniOption::IsotropicStress(value) => {
                state.stress.sym_set(0, 0, value);
                state.stress.sym_set(1, 1, value);
                state.stress.sym_set(2, 2, value);
            }
            IniOption::Zero => (),
        };
        */
        Ok(())
    }
}

// for point in &config.mesh.points {
//     if let Some(eq) = equation_numbers.get_option_equation_number(point.id, Dof::Pl) {
//         sim_state.system_xx[eq] = geostatics.calc_liquid_pressure(&point.coords)?;
//     }
// }
