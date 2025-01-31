use crate::material::{LocalState, LocalStatePorousLiq, LocalStatePorousSldLiq};
use crate::StrError;
use russell_tensor::{Mandel, Tensor2};
use serde::{Deserialize, Serialize};

/// Holds the secondary values (e.g., stress) at all integration points
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    /// Holds the number of integration points
    pub ngauss: usize,

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
            solid: Vec::new(),
            porous_liq: Vec::new(),
            porous_liq_gas: Vec::new(),
            porous_sld_liq: Vec::new(),
            porous_sld_liq_gas: Vec::new(),
        }
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
        if self.solid.len() == self.ngauss {
            Ok(&self.solid[p].stress)
        } else if self.porous_liq.len() == self.ngauss {
            Err("stress is not available in PorousLiq")
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
    /// config.update_model_settings(cell_attribute).save_strain = true;
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
        if self.solid.len() == self.ngauss {
            Ok(self.solid[p]
                .strain
                .as_ref()
                .ok_or("the recording of strains must be enabled first")?)
        } else if self.porous_liq.len() == self.ngauss {
            Err("strain is not available in PorousLiq")
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
}
