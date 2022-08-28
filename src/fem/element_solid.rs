#![allow(unused)]

use super::ElementEquations;
use crate::base::{Config, DofNumbers, ParamSolid};
use crate::StrError;
use gemlab::mesh::{Cell, Mesh};
use russell_tensor::Tensor2;

pub struct ElementSolid {
    pub sigma: Vec<Tensor2>, // (nip)
}

impl ElementSolid {
    pub fn new(
        mesh: &Mesh,
        dn: &DofNumbers,
        config: &Config,
        cell: &Cell,
        param: &ParamSolid,
    ) -> Result<Self, StrError> {
        let zero_sigma = Tensor2::new(true, true);
        Ok(ElementSolid {
            sigma: vec![zero_sigma; 1],
        })
    }
}

impl ElementEquations for ElementSolid {
    fn residual(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
    fn jacobian(&mut self) -> Result<(), StrError> {
        Err("stop")
    }
}
