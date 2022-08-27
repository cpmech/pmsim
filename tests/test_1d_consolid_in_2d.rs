#![allow(unused_variables)]
#![allow(unused_imports)]

use gemlab::mesh::{At, Extract, Features, Find};
use pmsim::base::{BcsEssential, BcsNatural, Dof, SampleMeshes};
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_quads();
    let features = Features::new(&mesh, Extract::Boundary);

    let find = Find::new(&mesh, &features)?;
    let left = find.edges(At::X(0.0))?;
    let right = find.edges(At::X(0.5))?;
    let bottom = find.edges(At::Y(0.0))?;
    let top = find.edges(At::Y(3.0))?;

    let zero = |_| 0.0;
    let mut bcs_essential = BcsEssential::new();
    bcs_essential
        .on(&left, &[Dof::Ux], zero)
        .on(&right, &[Dof::Ux], zero)
        .on(&bottom, &[Dof::Uy], zero);

    Ok(())
}
