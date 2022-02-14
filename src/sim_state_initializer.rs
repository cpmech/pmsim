use crate::{ModelRealDensity, ParamRealDensity, SimConfig, StateStress, StrError};
use russell_lab::Vector;

// /// Holds initialization data
// pub struct InitializationData {
//     /// At-rest earth pressure coefficient K0 = σₕ'/σᵥ' to compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
//     kk0: Option<f64>,

//     /// Poisson's coefficient ν to estimate the at-rest earth pressure coefficient K0 = ν/(1-ν) = σₕ'/σᵥ' and then compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
//     nu: Option<f64>,
// }

pub enum IniOption {
    /// Geostatic initial state where the parameter is the density of water
    ///
    /// # Note
    ///
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max (2D) or z=z_max (3D),
    ///   thus only fully water-saturated states are considered
    Geostatic(ParamRealDensity),

    /// Initial isotropic stress state where the parameter is σ_iso = σ_xx = σ_yy = σ_zz
    IsotropicStress(f64),

    /// Zero initial state
    Zero,
}

pub struct SimStateInitializer<'a> {
    config: &'a SimConfig<'a>, // Access to configuration
    model_water: Option<ModelRealDensity>,
}

impl<'a> SimStateInitializer<'a> {
    pub fn new(config: &'a SimConfig<'a>) -> Result<Self, StrError> {
        let model_water = match &config.ini_option {
            IniOption::Geostatic(params) => Some(ModelRealDensity::new(params)?),
            _ => None,
        };
        Ok(SimStateInitializer { config, model_water })
    }

    pub fn initialize_stress(&self, state: &mut StateStress, ip_coords: &Vector) -> Result<(), StrError> {
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
        Ok(())
    }
}

// for point in &config.mesh.points {
//     if let Some(eq) = equation_numbers.get_option_equation_number(point.id, Dof::Pl) {
//         sim_state.system_xx[eq] = geostatics.calc_liquid_pressure(&point.coords)?;
//     }
// }
