#![allow(unused)]

use crate::StrError;
use russell_lab::{Matrix, Vector};

pub struct Rod {}

impl Rod {
    pub fn fn_residual(residual: &mut Vector, time: f64, thickness: f64) -> Result<(), StrError> {
        Err("stop")
    }

    pub fn fn_jacobian(jacobian: &mut Matrix, time: f64, thickness: f64) -> Result<(), StrError> {
        Err("stop")
    }
}
