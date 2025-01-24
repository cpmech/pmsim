#![allow(unused)]

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_qua4();

    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(0.5), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(3.0), any_x)?;

    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .edges(&left, Ebc::Ux(zero))
        .edges(&right, Ebc::Ux(zero))
        .edges(&bottom, Ebc::Uy(zero));

    Ok(())
}
