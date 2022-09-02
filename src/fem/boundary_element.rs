#![allow(unused)]

use super::State;
use crate::base::{BcsNatural, DofNumbers, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Feature, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};

pub struct BoundaryElement {
    pub nbc: Nbc,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Option<Matrix>,
    pub local_to_global: Vec<usize>,
}
