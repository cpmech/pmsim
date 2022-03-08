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
    pub stress: Tensor2,           // (effective) stress
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

impl StateSeepage {
    /// Allocates a new instance
    pub fn new() -> Self {
        StateSeepage {
            ns0: 0.0,
            sat_liq: 0.0,
            real_rho_liq: 0.0,
            real_rho_gas: 0.0,
            delta_pc: 0.0,
        }
    }
}

impl StateStress {
    /// Allocates a new instance
    pub fn new(n_internal_values: usize, two_dim: bool) -> Self {
        StateStress {
            stress: Tensor2::new(true, two_dim),
            internal_values: vec![0.0; n_internal_values],
        }
    }
}

impl StateElement {
    /// Allocates a new instance with empty arrays
    pub fn new_empty() -> Self {
        StateElement {
            seepage: Vec::new(),
            stress: Vec::new(),
        }
    }

    /// Allocates a new instance with non-empty seepage array and with empty stress array
    pub fn new_seepage_only(n_integ_point: usize) -> Self {
        StateElement {
            seepage: vec![StateSeepage::new(); n_integ_point],
            stress: Vec::new(),
        }
    }

    /// Allocates a new instance with non-empty stress array and empty seepage array
    pub fn new_stress_only(n_integ_point: usize, n_internal_values: usize, two_dim: bool) -> Self {
        StateElement {
            seepage: Vec::new(),
            stress: vec![StateStress::new(n_internal_values, two_dim); n_integ_point],
        }
    }

    /// Allocates a new instance with non-empty seepage and stress arrays
    pub fn new_seepage_and_stress(n_integ_point: usize, n_internal_values: usize, two_dim: bool) -> Self {
        StateElement {
            seepage: vec![StateSeepage::new(); n_integ_point],
            stress: vec![StateStress::new(n_internal_values, two_dim); n_integ_point],
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
