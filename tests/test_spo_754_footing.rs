use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

const SAVE_FIGURE: bool = true;

#[test]
fn test_solid_smith_5d2_tri3_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::from_text_file("data/meshes/test_spo_754_footing.msh")?;
    if SAVE_FIGURE {
        mesh.check_all()?;
        let mut opt = Figure::new();
        opt.point_dots = true;
        // opt.point_ids = true;
        // opt.cell_ids = true;
        // opt.figure_size = Some((2000.0, 2000.0));
        opt.figure_size = Some((800.0, 800.0));
        mesh.draw(Some(opt), "/tmp/pmsim/test_spo_754_footing.svg", |_, _| {})?;
    }

    Ok(())
}
