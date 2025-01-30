use gemlab::mesh::Mesh;
use pmsim::fem::{FemBase, FemState, FileIo};
use pmsim::StrError;
use structopt::StructOpt;

/// Command line options
#[derive(StructOpt, Debug)]
#[structopt(
    name = "pmsim_to_paraview",
    about = "Generates VTU and PVD files for visualization with Paraview"
)]
struct Options {
    path_summary: String,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load FileIo
    let file_io = FileIo::read_json(&options.path_summary)?;

    // load the mesh
    let path_mesh = file_io.path_mesh();
    let mesh = Mesh::read_json(&path_mesh)?;

    // load the FemBase
    let path_base = file_io.path_base();
    let base = FemBase::read_json(&path_base)?;

    // write VTU files
    for index in &file_io.indices {
        // load state
        let path_state = file_io.path_state(*index);
        let state = FemState::read_json(&path_state)?;

        // write file
        file_io.write_vtu(&mesh, &base, &state, *index)?;
    }

    // write PVD file
    file_io.write_pvd()?;

    // message
    let path_pvd = file_io.path_pvd();
    let thin_line = format!("{:─^1$}", "", path_pvd.len());
    println!("\n\n{}", thin_line);
    println!("VTU files generated; the PVD file is:");
    println!("{}", path_pvd);
    println!("{}\n\n", thin_line);
    Ok(())
}
