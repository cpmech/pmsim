use pmsim::fem::PostProc;
use pmsim::StrError;
use structopt::StructOpt;

/// Command line options
#[derive(StructOpt)]
#[structopt(
    name = "pmsim2pv2d",
    about = "Generates VTU and PVD files for visualization with Paraview"
)]
struct Options {
    out_dir: String,

    fn_stem: String,

    #[structopt(long)]
    with_elastic_flags: bool,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load data
    let (post, mut memo) = PostProc::<2>::new(&options.out_dir, &options.fn_stem)?;

    // write VTU files
    for index in 0..post.nfile() {
        let state = post.read_file(index)?;
        post.write_vtu(
            &mut memo,
            &options.out_dir,
            &options.fn_stem,
            &state,
            index,
            options.with_elastic_flags,
        )?;
    }

    // write PVD file
    let path_pvd = post.write_pvd(&options.out_dir, &options.fn_stem)?;

    // message
    let thin_line = format!("{:─^1$}", "", path_pvd.len());
    println!("\n\n{}", thin_line);
    println!("VTU files generated; the PVD file is:");
    println!("{}", path_pvd);
    println!("{}\n\n", thin_line);
    Ok(())
}
