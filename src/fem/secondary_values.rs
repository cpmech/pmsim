use crate::material::{LocalState, LocalStatePorousLiq, LocalStatePorousSldLiq};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds the secondary values (e.g., stress) at all integration points
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    /// Holds the number of integration points
    pub ngauss: usize,

    /// Holds the flux vector at all integration points of a diffusion element
    ///
    /// **Note:** This field is used in post-processing only and must be enabled via [crate::base::Config::set_out_flux()].
    ///
    /// (ngauss)
    pub(crate) diffusion: Vec<Vector>,

    /// Holds the local states at all integration points of a solid element
    ///
    /// (ngauss)
    pub(crate) solid: Vec<LocalState>,

    /// Holds the local states at all integration points of a porous-liq element
    ///
    /// (ngauss)
    pub(crate) porous_liq: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-liq-gas element
    ///
    /// (ngauss)
    pub(crate) porous_liq_gas: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq element
    ///
    /// (ngauss)
    pub(crate) porous_sld_liq: Vec<LocalStatePorousSldLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq-gas element
    ///
    /// (ngauss)
    pub(crate) porous_sld_liq_gas: Vec<LocalStatePorousSldLiq>,
}

impl SecondaryValues {
    /// Allocates a new instance with empty arrays
    pub(crate) fn new_empty() -> Self {
        SecondaryValues {
            ngauss: 0,
            diffusion: Vec::new(),
            solid: Vec::new(),
            porous_liq: Vec::new(),
            porous_liq_gas: Vec::new(),
            porous_sld_liq: Vec::new(),
            porous_sld_liq_gas: Vec::new(),
        }
    }

    /// Allocates secondary values used in post-processing only of Diffusion elements
    pub(crate) fn allocate_diffusion(&mut self, ngauss: usize, ndim: usize) {
        let zero = Vector::new(ndim);
        self.diffusion = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Allocates secondary values for Solid elements
    pub(crate) fn allocate_solid(&mut self, mandel: Mandel, ngauss: usize, n_int_var: usize) {
        let zero = LocalState::new(mandel, n_int_var);
        self.solid = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Allocates secondary values for PorousLiq elements
    pub(crate) fn allocate_porous_liq(&mut self, ngauss: usize) {
        let zero = LocalStatePorousLiq::new();
        self.porous_liq = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Allocates secondary values for PorousLiqGas elements
    pub(crate) fn allocate_porous_liq_gas(&mut self, ngauss: usize) {
        let zero = LocalStatePorousLiq::new();
        self.porous_liq_gas = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Allocates secondary values for PorousSldLiq elements
    pub(crate) fn allocate_porous_sld_liq(&mut self, mandel: Mandel, ngauss: usize, n_int_var: usize) {
        let zero = LocalStatePorousSldLiq::new(mandel, n_int_var);
        self.porous_sld_liq = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Allocates secondary values for PorousSldLiqGas elements
    pub(crate) fn allocate_porous_sld_liq_gas(&mut self, mandel: Mandel, ngauss: usize, n_int_var: usize) {
        let zero = LocalStatePorousSldLiq::new(mandel, n_int_var);
        self.porous_sld_liq_gas = vec![zero; ngauss];
        self.ngauss = ngauss;
    }

    /// Returns the LocalState at an integration point
    ///
    /// # Input
    ///
    /// * `p` -- index of the integration point
    pub fn get_local_state(&self, p: usize) -> Result<&LocalState, StrError> {
        if self.ngauss == 0 {
            return Err("secondary values have not been allocated yet");
        }
        if p >= self.ngauss {
            return Err("index of integration point is out of bounds");
        }
        if self.diffusion.len() == self.ngauss {
            Err("LocalState is not available for Diffusion")
        } else if self.solid.len() == self.ngauss {
            Ok(&self.solid[p])
        } else if self.porous_liq.len() == self.ngauss {
            Err("LocalState is not available for PorousLiq")
        } else if self.porous_liq_gas.len() == self.ngauss {
            Err("LocalState is not available for PorousLiqGas")
        } else if self.porous_sld_liq.len() == self.ngauss {
            Err("LocalState is not available for PorousSldLiq")
        } else {
            Err("LocalState is not available")
        }
    }

    /// Returns the flux vector at an integration point
    ///
    /// # Input
    ///
    /// * `p` -- index of the integration point
    pub fn get_flux_vector(&self, p: usize) -> Result<&Vector, StrError> {
        if self.ngauss == 0 {
            return Err("secondary values have not been allocated yet");
        }
        if p >= self.ngauss {
            return Err("index of integration point is out of bounds");
        }
        if self.diffusion.len() == self.ngauss {
            Ok(&self.diffusion[p])
        } else if self.solid.len() == self.ngauss {
            Err("flow vector is not available for Solid")
        } else if self.porous_liq.len() == self.ngauss {
            Err("flow vector is not available for PorousLiq")
        } else if self.porous_liq_gas.len() == self.ngauss {
            Err("flow vector is not available for PorousLiqGas")
        } else if self.porous_sld_liq.len() == self.ngauss {
            Err("flow vector is not available for PorousSldLiq")
        } else {
            Err("flow vector is not available")
        }
    }

    /// Returns the stress tensor at an integration point
    ///
    /// # Input
    ///
    /// * `p` -- index of the integration point
    pub fn stress(&self, p: usize) -> Result<&Tensor2, StrError> {
        if self.ngauss == 0 {
            return Err("secondary values have not been allocated yet");
        }
        if p >= self.ngauss {
            return Err("index of integration point is out of bounds");
        }
        if self.diffusion.len() == self.ngauss {
            Err("stress is not available for Diffusion")
        } else if self.solid.len() == self.ngauss {
            Ok(&self.solid[p].stress)
        } else if self.porous_liq.len() == self.ngauss {
            Err("stress is not available for PorousLiq")
        } else if self.porous_liq_gas.len() == self.ngauss {
            Err("stress is not available for PorousLiqGas")
        } else if self.porous_sld_liq.len() == self.ngauss {
            Ok(&self.porous_sld_liq[p].stress)
        } else {
            Ok(&self.porous_sld_liq_gas[p].stress)
        }
    }

    /// Returns the strain tensor at an integration point
    ///
    /// Note: the recording of strains must be enabled in [crate::base::Config] first.
    /// For example:
    ///
    /// ````text
    /// config.update_model_settings(cell_marker).save_strain = true;
    /// ```
    ///
    /// # Input
    ///
    /// * `p` -- index of the integration point
    pub fn strain(&self, p: usize) -> Result<&Tensor2, StrError> {
        if self.ngauss == 0 {
            return Err("secondary values have not been allocated yet");
        }
        if p >= self.ngauss {
            return Err("index of integration point is out of bounds");
        }
        if self.diffusion.len() == self.ngauss {
            Err("strain is not available for Diffusion")
        } else if self.solid.len() == self.ngauss {
            Ok(self.solid[p]
                .strain
                .as_ref()
                .ok_or("the recording of strains must be enabled first")?)
        } else if self.porous_liq.len() == self.ngauss {
            Err("strain is not available for PorousLiq")
        } else if self.porous_liq_gas.len() == self.ngauss {
            Err("strain is not available for PorousLiqGas")
        } else if self.porous_sld_liq.len() == self.ngauss {
            Ok(self.porous_sld_liq[p]
                .strain
                .as_ref()
                .ok_or("the recording of strains must be enabled first")?)
        } else {
            Ok(self.porous_sld_liq_gas[p]
                .strain
                .as_ref()
                .ok_or("the recording of strains must be enabled first")?)
        }
    }

    /// Returns the elastic flag at an integration point
    ///
    /// # Input
    ///
    /// * `p` -- index of the integration point
    pub fn elastic_flag(&self, p: usize) -> Result<bool, StrError> {
        if self.ngauss == 0 {
            return Err("secondary values have not been allocated yet");
        }
        if p >= self.ngauss {
            return Err("index of integration point is out of bounds");
        }
        if self.diffusion.len() == self.ngauss {
            Err("elastic flag is not available for Diffusion")
        } else if self.solid.len() == self.ngauss {
            Ok(self.solid[p].elastic)
        } else if self.porous_liq.len() == self.ngauss {
            Err("elastic flag is not available for PorousLiq")
        } else if self.porous_liq_gas.len() == self.ngauss {
            Err("elastic flag is not available for PorousLiqGas")
        } else if self.porous_sld_liq.len() == self.ngauss {
            Ok(self.porous_sld_liq[p].elastic)
        } else {
            Ok(self.porous_sld_liq_gas[p].elastic)
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SecondaryValues;
    use russell_tensor::Mandel;

    #[test]
    fn new_empty_works() {
        let sv = SecondaryValues::new_empty();
        assert_eq!(sv.ngauss, 0);
        assert_eq!(sv.diffusion.len(), 0);
        assert_eq!(sv.solid.len(), 0);
        assert_eq!(sv.porous_liq.len(), 0);
        assert_eq!(sv.porous_liq_gas.len(), 0);
        assert_eq!(sv.porous_sld_liq.len(), 0);
        assert_eq!(sv.porous_sld_liq_gas.len(), 0);
    }

    #[test]
    fn allocate_diffusion_works() {
        let mut sv = SecondaryValues::new_empty();
        let ngauss = 4;
        let ndim = 2;

        sv.allocate_diffusion(ngauss, ndim);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.diffusion.len(), ngauss);
        for i in 0..ngauss {
            assert_eq!(sv.diffusion[i].dim(), ndim);
        }
    }

    #[test]
    fn allocate_solid_works() {
        let mut sv = SecondaryValues::new_empty();
        let mandel = Mandel::new(2);
        let ngauss = 4;
        let n_int_var = 2;

        sv.allocate_solid(mandel, ngauss, n_int_var);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.solid.len(), ngauss);
        for i in 0..ngauss {
            assert_eq!(sv.solid[i].int_vars.dim(), n_int_var);
        }
    }

    #[test]
    fn allocate_porous_liq_works() {
        let mut sv = SecondaryValues::new_empty();
        let ngauss = 4;

        sv.allocate_porous_liq(ngauss);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.porous_liq.len(), ngauss);
    }

    #[test]
    fn allocate_porous_liq_gas_works() {
        let mut sv = SecondaryValues::new_empty();
        let ngauss = 4;

        sv.allocate_porous_liq_gas(ngauss);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.porous_liq_gas.len(), ngauss);
    }

    #[test]
    fn allocate_porous_sld_liq_works() {
        let mut sv = SecondaryValues::new_empty();
        let mandel = Mandel::new(2);
        let ngauss = 4;
        let n_int_var = 3;

        sv.allocate_porous_sld_liq(mandel, ngauss, n_int_var);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.porous_sld_liq.len(), ngauss);
        for i in 0..ngauss {
            assert_eq!(sv.porous_sld_liq[i].int_vars.dim(), n_int_var);
        }
    }

    #[test]
    fn allocate_porous_sld_liq_gas_works() {
        let mut sv = SecondaryValues::new_empty();
        let mandel = Mandel::new(2);
        let ngauss = 4;
        let n_int_var = 3;

        sv.allocate_porous_sld_liq_gas(mandel, ngauss, n_int_var);

        assert_eq!(sv.ngauss, ngauss);
        assert_eq!(sv.porous_sld_liq_gas.len(), ngauss);
        for i in 0..ngauss {
            assert_eq!(sv.porous_sld_liq_gas[i].int_vars.dim(), n_int_var);
        }
    }

    #[test]
    fn get_local_state_handles_errors() {
        let sv = SecondaryValues::new_empty();

        // Not allocated yet
        assert_eq!(
            sv.get_local_state(0).err(),
            Some("secondary values have not been allocated yet")
        );

        // Allocate solid and test valid access
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 4, 2);

        // Valid access
        assert!(sv.get_local_state(0).is_ok());
        assert!(sv.get_local_state(3).is_ok());

        // Out of bounds
        assert_eq!(
            sv.get_local_state(4).err(),
            Some("index of integration point is out of bounds")
        );
    }

    #[test]
    fn get_local_state_works_for_solid() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 2, 1);

        let local_state = sv.get_local_state(0).unwrap();
        assert_eq!(local_state.int_vars.dim(), 1);
    }

    #[test]
    fn get_local_state_fails_for_diffusion() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_diffusion(2, 2);

        assert_eq!(
            sv.get_local_state(0).err(),
            Some("LocalState is not available for Diffusion")
        );
    }

    #[test]
    fn get_local_state_fails_for_porous_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq(2);

        assert_eq!(
            sv.get_local_state(0).err(),
            Some("LocalState is not available for PorousLiq")
        );
    }

    #[test]
    fn get_local_state_fails_for_porous_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq_gas(2);

        assert_eq!(
            sv.get_local_state(0).err(),
            Some("LocalState is not available for PorousLiqGas")
        );
    }

    #[test]
    fn get_local_state_fails_for_porous_sld_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq(Mandel::new(2), 2, 1);

        assert_eq!(
            sv.get_local_state(0).err(),
            Some("LocalState is not available for PorousSldLiq")
        );
    }

    #[test]
    fn get_local_state_fails_for_porous_sld_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq_gas(Mandel::new(2), 2, 1);

        assert_eq!(sv.get_local_state(0).err(), Some("LocalState is not available"));
    }

    #[test]
    fn get_flow_vector_handles_errors() {
        let sv = SecondaryValues::new_empty();

        // Not allocated yet
        assert_eq!(
            sv.get_flux_vector(0).err(),
            Some("secondary values have not been allocated yet")
        );

        // Allocate diffusion and test valid access
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_diffusion(4, 2);

        // Valid access
        assert!(sv.get_flux_vector(0).is_ok());
        assert!(sv.get_flux_vector(3).is_ok());

        // Out of bounds
        assert_eq!(
            sv.get_flux_vector(4).err(),
            Some("index of integration point is out of bounds")
        );
    }

    #[test]
    fn get_flow_vector_works_for_diffusion() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_diffusion(2, 3);

        let flow = sv.get_flux_vector(0).unwrap();
        assert_eq!(flow.dim(), 3);
    }

    #[test]
    fn get_flow_vector_fails_for_solid() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 2, 1);

        assert_eq!(
            sv.get_flux_vector(0).err(),
            Some("flow vector is not available for Solid")
        );
    }

    #[test]
    fn get_flow_vector_fails_for_porous_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq(2);

        assert_eq!(
            sv.get_flux_vector(0).err(),
            Some("flow vector is not available for PorousLiq")
        );
    }

    #[test]
    fn get_flow_vector_fails_for_porous_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq_gas(2);

        assert_eq!(
            sv.get_flux_vector(0).err(),
            Some("flow vector is not available for PorousLiqGas")
        );
    }

    #[test]
    fn get_flow_vector_fails_for_porous_sld_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq(Mandel::new(2), 2, 1);

        assert_eq!(
            sv.get_flux_vector(0).err(),
            Some("flow vector is not available for PorousSldLiq")
        );
    }

    #[test]
    fn get_flow_vector_fails_for_porous_sld_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq_gas(Mandel::new(2), 2, 1);

        assert_eq!(sv.get_flux_vector(0).err(), Some("flow vector is not available"));
    }

    #[test]
    fn stress_handles_errors() {
        let sv = SecondaryValues::new_empty();

        // Not allocated yet
        assert_eq!(sv.stress(0).err(), Some("secondary values have not been allocated yet"));

        // Allocate solid and test valid access
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 4, 2);

        // Valid access
        assert!(sv.stress(0).is_ok());
        assert!(sv.stress(3).is_ok());

        // Out of bounds
        assert_eq!(sv.stress(4).err(), Some("index of integration point is out of bounds"));
    }

    #[test]
    fn stress_works_for_solid() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 2, 1);

        let stress = sv.stress(0).unwrap();
        assert_eq!(stress.mandel(), Mandel::new(2));
    }

    #[test]
    fn stress_works_for_porous_sld_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq(Mandel::new(2), 2, 1);

        let stress = sv.stress(0).unwrap();
        assert_eq!(stress.mandel(), Mandel::new(2));
    }

    #[test]
    fn stress_works_for_porous_sld_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq_gas(Mandel::new(2), 2, 1);

        let stress = sv.stress(0).unwrap();
        assert_eq!(stress.mandel(), Mandel::new(2));
    }

    #[test]
    fn stress_fails_for_diffusion() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_diffusion(2, 2);

        assert_eq!(sv.stress(0).err(), Some("stress is not available for Diffusion"));
    }

    #[test]
    fn stress_fails_for_porous_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq(2);

        assert_eq!(sv.stress(0).err(), Some("stress is not available for PorousLiq"));
    }

    #[test]
    fn stress_fails_for_porous_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq_gas(2);

        assert_eq!(sv.stress(0).err(), Some("stress is not available for PorousLiqGas"));
    }

    #[test]
    fn strain_handles_errors() {
        let sv = SecondaryValues::new_empty();

        // Not allocated yet
        assert_eq!(sv.strain(0).err(), Some("secondary values have not been allocated yet"));

        // Allocate solid and test (without strain recording enabled)
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 4, 2);

        // Out of bounds
        assert_eq!(sv.strain(4).err(), Some("index of integration point is out of bounds"));

        // Strain not enabled (default)
        assert_eq!(
            sv.strain(0).err(),
            Some("the recording of strains must be enabled first")
        );
    }

    #[test]
    fn strain_works_for_solid_when_enabled() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_solid(Mandel::new(2), 2, 1);

        // Enable strain recording
        sv.solid[0].enable_strain();
        sv.solid[1].enable_strain();

        let strain = sv.strain(0).unwrap();
        assert_eq!(strain.mandel(), Mandel::new(2));
    }

    #[test]
    fn strain_works_for_porous_sld_liq_when_enabled() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq(Mandel::new(2), 2, 1);

        // Enable strain recording
        sv.porous_sld_liq[0].enable_strain();

        let strain = sv.strain(0).unwrap();
        assert_eq!(strain.mandel(), Mandel::new(2));
    }

    #[test]
    fn strain_works_for_porous_sld_liq_gas_when_enabled() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_sld_liq_gas(Mandel::new(2), 2, 1);

        // Enable strain recording
        sv.porous_sld_liq_gas[0].enable_strain();

        let strain = sv.strain(0).unwrap();
        assert_eq!(strain.mandel(), Mandel::new(2));
    }

    #[test]
    fn strain_fails_for_diffusion() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_diffusion(2, 2);

        assert_eq!(sv.strain(0).err(), Some("strain is not available for Diffusion"));
    }

    #[test]
    fn strain_fails_for_porous_liq() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq(2);

        assert_eq!(sv.strain(0).err(), Some("strain is not available for PorousLiq"));
    }

    #[test]
    fn strain_fails_for_porous_liq_gas() {
        let mut sv = SecondaryValues::new_empty();
        sv.allocate_porous_liq_gas(2);

        assert_eq!(sv.strain(0).err(), Some("strain is not available for PorousLiqGas"));
    }
}
