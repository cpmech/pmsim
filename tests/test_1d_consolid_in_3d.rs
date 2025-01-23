#![allow(unused)]

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_qua4();

    let features = Features::new(&mesh, false);
    let left = features.search_faces(At::X(0.0), any_x)?;
    let right = features.search_faces(At::X(0.5), any_x)?;
    let back = features.search_faces(At::X(0.0), any_x)?;
    let front = features.search_faces(At::X(0.5), any_x)?;
    let bottom = features.search_faces(At::Y(0.0), any_x)?;
    let top = features.search_faces(At::Y(3.0), any_x)?;

    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(zero))
        .on(&right, Ebc::Ux(zero))
        .on(&bottom, Ebc::Uy(zero));

    Ok(())
}
