use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct StateSeepage {
    ns0: f64,          // initial partial fraction of solids
    sat_liq: f64,      // liquid saturation
    real_rho_liq: f64, // real (intrinsic) density of liquid
    real_rho_gas: f64, // real (intrinsic) density of gas
    delta_pc: f64,     // step increment of capillary pressure
    wetting: bool,     // wetting flag
}
