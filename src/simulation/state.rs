use russell_lab::Vector;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

/// Holds seepage state variables for an integration point
#[derive(Clone, Debug, Deserialize, Serialize)]
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

/// Holds all simulation state variables
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct State {
    pub elements: Vec<StateElement>, // (nele)
    pub system_xx: Vector,           // (neq)
    pub system_yy: Vector,           // (neq)
}

impl StateElement {
    /// Allocates a new instance with empty arrays
    pub fn new_empty() -> Self {
        StateElement {
            seepage: Vec::new(),
            stress: Vec::new(),
        }
    }
}

impl State {
    /// Allocate a new instance with empty arrays
    pub fn new_empty() -> Self {
        State {
            elements: Vec::new(),
            system_xx: Vector::new(0),
            system_yy: Vector::new(0),
        }
    }
}
