#![allow(unused)]

use super::ElementEquations;
use crate::base::ParamSolid;
use crate::StrError;
use russell_lab::{Matrix, Vector};
use russell_tensor::Tensor2;

pub struct Solid {
    pub sigma: Vec<Tensor2>, // (nip)
}

impl Solid {
    pub fn new(param: ParamSolid) -> Self {
        let zero_sigma = Tensor2::new(true, true);
        Solid {
            sigma: vec![zero_sigma; 1],
        }
    }

    pub fn fn_residual(residual: &mut Vector, time: f64, thickness: f64) -> Result<(), StrError> {
        Err("stop")
    }

    pub fn fn_jacobian(jacobian: &mut Matrix, time: f64, thickness: f64) -> Result<(), StrError> {
        Err("stop")
    }
}

impl ElementEquations for Solid {
    fn residual(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
    fn jacobian(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
}
