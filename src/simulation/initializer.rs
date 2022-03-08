use super::Configuration;
use crate::{geostatics::Geostatics, StrError};
use russell_lab::Vector;
use russell_tensor::Tensor2;

pub enum IniOption {
    /// Geostatic initial state with data = (overburden,total_stress)
    ///
    /// * The argument is the overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial isotropic stress state with σ_xx = σ_yy = σ_zz = value
    Isotropic(f64),

    /// Zero initial state
    Zero,
}

pub struct Initializer {
    space_ndim: usize,              // number of space dimensions
    total_stress: bool,             // total stress analysis?
    geostatics: Option<Geostatics>, // geostatic initialization
    isotropic: f64,                 // isotropic initialization
}

impl Initializer {
    /// Allocates a new instance
    pub fn new(config: &Configuration) -> Result<Self, StrError> {
        let space_ndim = config.get_mesh().space_ndim;
        let total_stress = config.get_total_stress();
        match config.get_ini_option() {
            &IniOption::Geostatic(..) => Ok(Initializer {
                space_ndim,
                total_stress,
                geostatics: Some(Geostatics::new(config)?),
                isotropic: 0.0,
            }),
            &IniOption::Isotropic(value) => Ok(Initializer {
                space_ndim,
                total_stress,
                geostatics: None,
                isotropic: value,
            }),
            &IniOption::Zero => Ok(Initializer {
                space_ndim,
                total_stress,
                geostatics: None,
                isotropic: 0.0,
            }),
        }
    }

    /// Returns total or effective stress at an integration point
    pub fn stress_at_ip(&self, ip_coords: &Vector) -> Result<Tensor2, StrError> {
        if ip_coords.dim() != self.space_ndim {
            return Err("ip_coords.dim() must be equal to space_ndim");
        }
        if let Some(geostatics) = &self.geostatics {
            let z = ip_coords[self.space_ndim - 1];
            return geostatics.calc_stress(z, self.total_stress);
        }
        let mut sigma = Tensor2::new(true, self.space_ndim == 2);
        sigma.vec[0] = self.isotropic;
        sigma.vec[1] = self.isotropic;
        sigma.vec[2] = self.isotropic;
        Ok(sigma)
    }
}

// for point in &config.mesh.points {
//     if let Some(eq) = equation_numbers.get_option_equation_number(point.id, Dof::Pl) {
//         sim_state.system_xx[eq] = geostatics.calc_liquid_pressure(&point.coords)?;
//     }
// }
