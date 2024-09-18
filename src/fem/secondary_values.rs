use crate::base::Idealization;
use crate::material::{LocalState, LocalStatePorousLiq, LocalStatePorousSldLiq};
use russell_tensor::Mandel;
use serde::{Deserialize, Serialize};

/// Holds the secondary values (e.g., stress) at all integration points
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct SecondaryValues {
    /// Holds the local states at all integration points of a solid element
    ///
    /// (n_integration_point)
    pub solid: Vec<LocalState>,

    /// Holds the local states at all integration points of a porous-liq element
    ///
    /// (n_integration_point)
    pub porous_liq: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-liq-gas element
    ///
    /// (n_integration_point)
    pub porous_liq_gas: Vec<LocalStatePorousLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq element
    ///
    /// (n_integration_point)
    pub porous_sld_liq: Vec<LocalStatePorousSldLiq>,

    /// Holds the local states at all integration points of a porous-sld-liq-gas element
    ///
    /// (n_integration_point)
    pub porous_sld_liq_gas: Vec<LocalStatePorousSldLiq>,

    /// (backup) Holds the local states at all integration points of a solid element
    bkp_solid: Vec<LocalState>,

    /// (backup) Holds the local states at all integration points of a porous-liq element
    bkp_porous_liq: Vec<LocalStatePorousLiq>,

    /// (backup) Holds the local states at all integration points of a porous-liq-gas element
    bkp_porous_liq_gas: Vec<LocalStatePorousLiq>,

    /// (backup) Holds the local states at all integration points of a porous-sld-liq element
    bkp_porous_sld_liq: Vec<LocalStatePorousSldLiq>,

    /// (backup) Holds the local states at all integration points of a porous-sld-liq-gas element
    bkp_porous_sld_liq_gas: Vec<LocalStatePorousSldLiq>,

    /// Holds the Mandel representation
    mandel: Mandel,
}

impl SecondaryValues {
    /// Allocates a new instance with empty arrays
    pub(crate) fn new_empty(ideal: &Idealization) -> Self {
        SecondaryValues {
            solid: Vec::new(),
            porous_liq: Vec::new(),
            porous_liq_gas: Vec::new(),
            porous_sld_liq: Vec::new(),
            porous_sld_liq_gas: Vec::new(),
            bkp_solid: Vec::new(),
            bkp_porous_liq: Vec::new(),
            bkp_porous_liq_gas: Vec::new(),
            bkp_porous_sld_liq: Vec::new(),
            bkp_porous_sld_liq_gas: Vec::new(),
            mandel: ideal.mandel(),
        }
    }

    /// Allocates secondary values for Solid elements
    pub(crate) fn allocate_solid(&mut self, n_integration_point: usize, n_internal_values: usize) {
        let zero = LocalState::new(self.mandel, n_internal_values);
        let bkp = zero.clone();
        self.solid = vec![zero; n_integration_point];
        self.bkp_solid = vec![bkp; n_integration_point];
    }

    /// Allocates secondary values for PorousLiq elements
    pub(crate) fn allocate_porous_liq(&mut self, n_integration_point: usize) {
        let zero = LocalStatePorousLiq::new();
        let bkp = zero.clone();
        self.porous_liq = vec![zero; n_integration_point];
        self.bkp_porous_liq = vec![bkp; n_integration_point];
    }

    /// Allocates secondary values for PorousLiqGas elements
    pub(crate) fn allocate_porous_liq_gas(&mut self, n_integration_point: usize) {
        let zero = LocalStatePorousLiq::new();
        let bkp = zero.clone();
        self.porous_liq_gas = vec![zero; n_integration_point];
        self.bkp_porous_liq_gas = vec![bkp; n_integration_point];
    }

    /// Allocates secondary values for PorousSldLiq elements
    pub(crate) fn allocate_porous_sld_liq(&mut self, n_integration_point: usize, n_internal_values: usize) {
        let zero = LocalStatePorousSldLiq::new(self.mandel, n_internal_values);
        let bkp = zero.clone();
        self.porous_sld_liq = vec![zero; n_integration_point];
        self.bkp_porous_sld_liq = vec![bkp; n_integration_point];
    }

    /// Allocates secondary values for PorousSldLiqGas elements
    pub(crate) fn allocate_porous_sld_liq_gas(&mut self, n_integration_point: usize, n_internal_values: usize) {
        let zero = LocalStatePorousSldLiq::new(self.mandel, n_internal_values);
        let bkp = zero.clone();
        self.porous_sld_liq_gas = vec![zero; n_integration_point];
        self.bkp_porous_sld_liq_gas = vec![bkp; n_integration_point];
    }

    /// Creates a copy of the secondary values (e.g., stresses and internal values)
    pub(crate) fn backup(&mut self) {
        self.bkp_solid
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.solid[i]));
        self.bkp_porous_liq
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.porous_liq[i]));
        self.bkp_porous_liq_gas
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.porous_liq_gas[i]));
        self.bkp_porous_sld_liq
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.porous_sld_liq[i]));
        self.bkp_porous_sld_liq_gas
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.porous_sld_liq_gas[i]));
    }

    /// Restores the secondary values from the backup (e.g., stresses and internal values)
    pub(crate) fn restore(&mut self) {
        self.solid
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.bkp_solid[i]));
        self.porous_liq
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.bkp_porous_liq[i]));
        self.porous_liq_gas
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.bkp_porous_liq_gas[i]));
        self.porous_sld_liq
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.bkp_porous_sld_liq[i]));
        self.porous_sld_liq_gas
            .iter_mut()
            .enumerate()
            .for_each(|(i, s)| s.mirror(&self.bkp_porous_sld_liq_gas[i]));
    }

    /// Resets algorithmic variables such as Λ at the beginning of implicit iterations
    pub(crate) fn reset_algorithmic_variables(&mut self) {
        self.solid //
            .iter_mut()
            .for_each(|s| s.reset_algorithmic_variables());
        self.porous_liq //
            .iter_mut()
            .for_each(|s| s.reset_algorithmic_variables());
        self.porous_liq_gas
            .iter_mut()
            .for_each(|s| s.reset_algorithmic_variables());
        self.porous_sld_liq
            .iter_mut()
            .for_each(|s| s.reset_algorithmic_variables());
        self.porous_sld_liq_gas
            .iter_mut()
            .for_each(|s| s.reset_algorithmic_variables());
    }
}
