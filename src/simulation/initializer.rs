use super::{Configuration, Dof, EquationId, State, UNASSIGNED};
use crate::{geostatics::Geostatics, StrError};
use gemlab::mesh::Mesh;
use russell_tensor::Tensor2;

/// Holds an option to initialize stresses
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

/// Implements functions to initialize stresses
pub struct Initializer {
    space_ndim: usize,              // number of space dimensions
    total_stress: bool,             // total stress analysis?
    geostatics: Option<Geostatics>, // geostatic initialization
    isotropic: f64,                 // isotropic initialization
}

impl Initializer {
    /// Allocates a new instance
    pub fn new(config: &Configuration) -> Result<Self, StrError> {
        let space_ndim = config.mesh.space_ndim;
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

    /// Returns the liquid pressure at given coordinates
    pub fn pl(&self, coords: &[f64]) -> Result<f64, StrError> {
        if coords.len() != self.space_ndim {
            return Err("coords.len() must be equal to space_ndim to calculate pl");
        }
        if let Some(geostatics) = &self.geostatics {
            let z = coords[self.space_ndim - 1];
            return geostatics.calc_pl(z);
        }
        Ok(0.0)
    }

    /// Returns the gas pressure at given coordinates
    pub fn pg(&self, coords: &[f64]) -> Result<f64, StrError> {
        if coords.len() != self.space_ndim {
            return Err("coords.len() must be equal to space_ndim to calculate pg");
        }
        if let Some(geostatics) = &self.geostatics {
            let z = coords[self.space_ndim - 1];
            return geostatics.calc_pg(z);
        }
        Ok(0.0)
    }

    /// Returns total or effective stress at given coordinates
    pub fn stress(&self, coords: &[f64]) -> Result<Tensor2, StrError> {
        if coords.len() != self.space_ndim {
            return Err("coords.len() must be equal to space_ndim to calculate stress");
        }
        if let Some(geostatics) = &self.geostatics {
            let z = coords[self.space_ndim - 1];
            return geostatics.calc_stress(z, self.total_stress);
        }
        let mut sigma = Tensor2::new(true, self.space_ndim == 2);
        sigma.vec[0] = self.isotropic;
        sigma.vec[1] = self.isotropic;
        sigma.vec[2] = self.isotropic;
        Ok(sigma)
    }

    /// Initializes liquid and/or gas pressure values
    pub fn liquid_gas_pressure(
        &self,
        state: &mut State,
        mesh: &Mesh,
        equation_id: &EquationId,
    ) -> Result<(), StrError> {
        for point in &mesh.points {
            // liquid pressure
            let eid = equation_id.eid(point.id, Dof::Pl);
            if eid != UNASSIGNED {
                let e = if eid < 0 { -eid } else { eid };
                let i = (e as usize) - 1;
                state.unknowns[i] = self.pl(&point.coords)?;
            }
            // gas pressure
            let eid = equation_id.eid(point.id, Dof::Pg);
            if eid != UNASSIGNED {
                let e = if eid < 0 { -eid } else { eid };
                let i = (e as usize) - 1;
                state.unknowns[i] = self.pg(&point.coords)?;
            }
        }
        Ok(())
    }
}
