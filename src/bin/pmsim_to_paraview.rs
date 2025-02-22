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
    let (post, _) = PostProc::new(&options.out_dir, &options.fn_stem)?;

    // write VTU files
    for index in 0..post.n_state() {
        let state = post.read_state(index)?;
        post.write_vtu(&state, index)?;
    }

    // write PVD file
    let path_pvd = post.write_pvd()?;

    // message
    let thin_line = format!("{:─^1$}", "", path_pvd.len());
    println!("\n\n{}", thin_line);
    println!("VTU files generated; the PVD file is:");
    println!("{}", path_pvd);
    println!("{}\n\n", thin_line);
    Ok(())
}
