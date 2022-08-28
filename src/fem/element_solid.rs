#![allow(unused)]

use super::ElementEquations;
use crate::base::{Config, DofNumbers, ParamSolid};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell, Mesh};
use gemlab::shapes::Scratchpad;
use russell_tensor::Tensor2;

pub struct ElementSolid {
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
}

impl ElementSolid {
    pub fn new(
        mesh: &Mesh,
        dn: &DofNumbers,
        config: &Config,
        cell: &Cell,
        param: &ParamSolid,
    ) -> Result<Self, StrError> {
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(mesh.ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, &mesh);
        Ok({
            ElementSolid {
                pad,
                ips: config.integ_point_data(cell)?,
            }
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
