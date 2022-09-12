#![allow(unused)]

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;

fn any(_: &Vec<f64>) -> bool {
    true
}

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_quads();
    let features = Features::new(&mesh, Extract::Boundary);

    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(0.0), any)?;
    let right = find.edges(At::X(0.5), any)?;
    let bottom = find.edges(At::Y(0.0), any)?;
    let top = find.edges(At::Y(3.0), any)?;

    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(zero))
        .on(&right, Ebc::Ux(zero))
        .on(&bottom, Ebc::Uy(zero));

    Ok(())
}
