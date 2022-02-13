use crate::{SimConfig, StateStress, StrError};
use russell_lab::Vector;

// /// Holds initialization data
// pub struct InitializationData {
//     /// At-rest earth pressure coefficient K0 = σₕ'/σᵥ' to compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
//     kk0: Option<f64>,

//     /// Poisson's coefficient ν to estimate the at-rest earth pressure coefficient K0 = ν/(1-ν) = σₕ'/σᵥ' and then compute horizontal effective stress (σₕ') from vertical effective stress (σᵥ')
//     nu: Option<f64>,
// }

pub enum IniOption {
    /// Geostatic initial state
    Geostatic,

    /// Initial isotropic stress state where the parameter is σ_iso = σ_xx = σ_yy = σ_zz
    IsotropicStress(f64),

    /// Zero initial state
    Zero,
}

pub struct SimStateInitializer<'a> {
    /// Access to configuration
    config: &'a SimConfig<'a>,
}

impl<'a> SimStateInitializer<'a> {
    pub fn new(config: &'a SimConfig<'a>) -> Result<Self, StrError> {
        // let initializer = match config.ini_option {
        //     IniOption::Geostatic => Geostatics::new(config)?,
        //     IniOption::IsotropicStress(..) => (),
        //     IniOption::Zero => (),
        // };

        Ok(SimStateInitializer { config })
    }

    pub fn calc_stress(&self, state: &mut StateStress, ip_coords: &Vector) -> Result<(), StrError> {
        Ok(())
    }
}

// for point in &config.mesh.points {
//     if let Some(eq) = equation_numbers.get_option_equation_number(point.id, Dof::Pl) {
//         sim_state.system_xx[eq] = geostatics.calc_liquid_pressure(&point.coords)?;
//     }
// }
