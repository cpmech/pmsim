#![allow(unused_variables)]
#![allow(unused_imports)]

use gemlab::mesh::{At, Extract, Features, Find};
use pmsim::base::{BcEssential, BcNatural, Dof, SampleMeshes};
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    let mesh = SampleMeshes::column_two_layers_quads();
    let features = Features::new(&mesh, Extract::Boundary);

    let find = Find::new(&mesh, &features)?;
    let left = find.faces(At::X(0.0))?;
    let right = find.faces(At::X(0.5))?;
    let back = find.faces(At::X(0.0))?;
    let front = find.faces(At::X(0.5))?;
    let bottom = find.faces(At::Y(0.0))?;
    let top = find.faces(At::Y(3.0))?;

    let zero = |_| 0.0;
    let mut ebc = BcEssential::new();
    ebc.set_faces(&left, &[Dof::Ux], zero)
        .set_faces(&right, &[Dof::Ux], zero)
        .set_faces(&bottom, &[Dof::Uy], zero);

    Ok(())
}
