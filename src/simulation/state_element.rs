use russell_tensor::{copy_tensor2, Tensor2};
use serde::{Deserialize, Serialize};

use crate::StrError;

/// Holds seepage state variables for an integration point
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub struct StateSeepage {
    pub ns0: f64,          // initial partial fraction of solids
    pub sat_liq: f64,      // liquid saturation sl
    pub real_rho_liq: f64, // real (intrinsic) density of liquid
    pub real_rho_gas: f64, // real (intrinsic) density of gas
    pub delta_pc: f64,     // step increment of capillary pressure
}

/// Holds stress state variables for an integration point
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct StateStress {
    pub sigma: Tensor2,            // total or effective stress
    pub internal_values: Vec<f64>, // internal values
}

/// Holds seepage and/or stress state variables for a set of integration points
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct StateElement {
    pub seepage: Vec<StateSeepage>, // (n_integ_point)
    pub stress: Vec<StateStress>,   // (n_integ_point)
}

impl StateElement {
    /// Allocates a new instance with empty arrays
    pub fn new_empty() -> Self {
        StateElement {
            seepage: Vec::new(),
            stress: Vec::new(),
        }
    }

    /// Copy all values into another StateElement
    pub fn copy_into(&self, other: &mut StateElement) -> Result<(), StrError> {
        if other.seepage.len() != self.seepage.len() {
            return Err("other StateElement has an incorrect number of seepage integration points");
        }
        if other.stress.len() != self.stress.len() {
            return Err("other StateElement has an incorrect number of stress integration points");
        }
        for i in 0..self.seepage.len() {
            other.seepage[i] = self.seepage[i];
        }
        for i in 0..self.stress.len() {
            copy_tensor2(&mut other.stress[i].sigma, &self.stress[i].sigma)?;
            other.stress[i]
                .internal_values
                .copy_from_slice(&self.stress[i].internal_values);
        }
        Ok(())
    }
}
