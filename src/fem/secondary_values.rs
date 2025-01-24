use crate::material::{LocalState, LocalStatePorousLiq, LocalStatePorousSldLiq};
use russell_tensor::Mandel;
use serde::{Deserialize, Serialize};

/// Holds the secondary values (e.g., stress) at all integration points
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    /// Holds the local states at all integration points of a solid element
    ///
    /// (ngauss)
    pub solid: Vec<LocalState>,

    /// Holds the local states at all integration points of a porous-liq element
    ///
    /// (ngauss)
    pub porous_liq: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-liq-gas element
    ///
    /// (ngauss)
    pub porous_liq_gas: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq element
    ///
    /// (ngauss)
    pub porous_sld_liq: Vec<LocalStatePorousSldLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq-gas element
    ///
    /// (ngauss)
    pub porous_sld_liq_gas: Vec<LocalStatePorousSldLiq>,
}

impl SecondaryValues {
    /// Allocates a new instance with empty arrays
    pub(crate) fn new_empty() -> Self {
        SecondaryValues {
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
    }

    /// Allocates secondary values for PorousLiq elements
    pub(crate) fn allocate_porous_liq(&mut self, ngauss: usize) {
        let zero = LocalStatePorousLiq::new();
        self.porous_liq = vec![zero; ngauss];
    }

    /// Allocates secondary values for PorousLiqGas elements
    pub(crate) fn allocate_porous_liq_gas(&mut self, ngauss: usize) {
        let zero = LocalStatePorousLiq::new();
        self.porous_liq_gas = vec![zero; ngauss];
    }

    /// Allocates secondary values for PorousSldLiq elements
    pub(crate) fn allocate_porous_sld_liq(&mut self, mandel: Mandel, ngauss: usize, n_int_var: usize) {
        let zero = LocalStatePorousSldLiq::new(mandel, n_int_var);
        self.porous_sld_liq = vec![zero; ngauss];
    }

    /// Allocates secondary values for PorousSldLiqGas elements
    pub(crate) fn allocate_porous_sld_liq_gas(&mut self, mandel: Mandel, ngauss: usize, n_int_var: usize) {
        let zero = LocalStatePorousSldLiq::new(mandel, n_int_var);
        self.porous_sld_liq_gas = vec![zero; ngauss];
    }
}
