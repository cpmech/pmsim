#![allow(unused)]

use crate::base::{BcEssential, BcNatural, Config, Schema};
use gemlab::mesh::Mesh;

pub struct Solver<'a> {
    mesh: &'a Mesh,
    schema: &'a Schema,
    config: &'a Config<'a>,
    essential: &'a BcEssential<'a>,
    // natural: &'a BcNatural<'a>,
}
