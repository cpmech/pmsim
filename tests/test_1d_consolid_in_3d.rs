#![allow(unused)]

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_qua4();

    let feat = Features::new(&mesh, false);
    let left = feat.search_faces(At::X(0.0), any_x)?;
    let right = feat.search_faces(At::X(0.5), any_x)?;
    let back = feat.search_faces(At::X(0.0), any_x)?;
    let front = feat.search_faces(At::X(0.5), any_x)?;
    let bottom = feat.search_faces(At::Y(0.0), any_x)?;
    let top = feat.search_faces(At::Y(3.0), any_x)?;

    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(zero))
        .on(&right, Ebc::Ux(zero))
        .on(&bottom, Ebc::Uy(zero));

    Ok(())
}
