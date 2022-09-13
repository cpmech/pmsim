#![allow(unused)]

use crate::base::ParamConductivity;
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::{copy_tensor2, Tensor2};

pub trait Conductivity: Send + Sync {
    /// Computes the conductivity tensor
    ///
    /// * `phi` may be the temperature or liquid/gas pressure
    fn tensor(&mut self, kk: &mut Tensor2, phi: f64, ivs: &Vector) -> Result<(), StrError>;
}

/// Allocates conductivity model
pub fn allocate_conductivity_model(param: &ParamConductivity, two_dim: bool) -> Box<dyn Conductivity> {
    let model: Box<dyn Conductivity> = match param {
        ParamConductivity::Constant { kx, ky, kz } => Box::new(ConductivityConstant::new(*kx, *ky, *kz, two_dim)),
        _ => panic!("todo"),
    };
    model
}

pub struct ConductivityConstant {
    pub kk: Tensor2,
}

impl ConductivityConstant {
    pub fn new(kx: f64, ky: f64, kz: f64, two_dim: bool) -> Self {
        let mut kk = Tensor2::new(true, two_dim);
        kk.sym_set(0, 0, kx);
        kk.sym_set(1, 1, ky);
        if !two_dim {
            kk.sym_set(2, 2, kz);
        }
        ConductivityConstant { kk }
    }
}

impl Conductivity for ConductivityConstant {
    fn tensor(&mut self, kk: &mut Tensor2, phi: f64, ivs: &Vector) -> Result<(), StrError> {
        copy_tensor2(kk, &self.kk)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::allocate_conductivity_model;

    #[test]
    fn allocate_conductivity_model_works() {}
}
