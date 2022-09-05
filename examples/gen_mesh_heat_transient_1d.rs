use gemlab::{prelude::*, StrError};

fn main() -> Result<(), StrError> {
    let mut block = Block::new(&[[0.0, 0.0], [20.0, 0.0], [20.0, 1.0], [0.0, 1.0]])?;
    block.set_ndiv(&[10, 1])?;
    let mesh = block.subdivide(GeoKind::Qua4)?;
    mesh.write("/tmp/pmsim/mesh_heat_transient_1d.dat")?;
    draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_transient_1d.svg")?;
    Ok(())
}
