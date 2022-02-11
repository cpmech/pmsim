#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::ElementConfig::*;
use crate::*;
use russell_lab::Vector;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct SimState {
    pub system_xx: Vector,                                    // (neq)
    pub system_yy: Vector,                                    // (neq)
    pub states_stress: Vec<Vec<StateStress>>,                 // (nele, nip_ele)
    pub states_seepage_liq: Vec<Vec<StateSeepageLiq>>,        // (nele, nip_ele)
    pub states_seepage_liq_gas: Vec<Vec<StateSeepageLiqGas>>, // (nele, nip_ele)
}

impl SimState {
    pub fn new(config: &SimConfig, neq: usize) -> Result<Self, StrError> {
        for cell in &config.mesh.cells {
            // get element configuration
            let element_config = match config.element_configs.get(&cell.attribute_id) {
                Some(c) => c,
                None => return Err("cannot find element configuration for a cell attribute id"),
            };
            // todo
            match element_config {
                Rod(..) => (),
                Beam(..) => (),
                Solid(params, nip) => {
                    let model = new_stress_strain_model(&params.stress_strain, config.two_dim, config.plane_stress);
                }
                SeepageLiq(params, nip) => (),
                SeepageLiqGas(params, nip) => (),
                PorousSolLiq(params, nip) => (),
                PorousSolLiqGas(params, nip) => (),
            }
        }
        Ok(SimState {
            system_xx: Vector::new(neq),
            system_yy: Vector::new(neq),
            states_stress: Vec::new(),
            states_seepage_liq: Vec::new(),
            states_seepage_liq_gas: Vec::new(),
        })
    }
}
