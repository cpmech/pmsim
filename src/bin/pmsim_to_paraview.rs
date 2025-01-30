use pmsim::fem::PostProc;
use pmsim::StrError;
use structopt::StructOpt;

/// Command line options
#[derive(StructOpt, Debug)]
#[structopt(
    name = "pmsim_to_paraview",
    about = "Generates VTU and PVD files for visualization with Paraview"
)]
struct Options {
    out_dir: String,

    fn_stem: String,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load data
    let (file_io, mesh, base) = PostProc::read_essential(&options.out_dir, &options.fn_stem)?;

    // write VTU files
    for index in &file_io.indices {
        let state = PostProc::read_state(&file_io, *index)?;
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
