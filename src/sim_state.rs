#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{StateStress, StrError};
use russell_lab::Vector;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SimState {
    pub system_xx: Vector,               // (neq)
    pub system_yy: Vector,               // (neq)
    pub stresses: Vec<Vec<StateStress>>, // (nele, nip_ele)
}

impl SimState {
    pub fn new(space_ndim: usize, neq: usize) -> Self {
        SimState {
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
            stresses: Vec::new(),
        }
    }

    pub fn initialize_stress() -> Result<(), StrError> {
        Ok(())
    }
}
