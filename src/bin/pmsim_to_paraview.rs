use pmsim::fem::FemState;
use pmsim::util::paraview_write_vtu;
use pmsim::StrError;
use structopt::StructOpt;

/// Command line options
#[derive(StructOpt, Debug)]
#[structopt(
    name = "pmsim_to_paraview",
    about = "Generates VTU files for visualization with Paraview"
)]
struct Options {
    /// Directory containing the results (.state) files
    results_directory: String,

    /// Filename stem
    filename_stem: String,

    /// Index of the output
    output_count: usize,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load state
    let path_state = format!(
        "{}/{}-{:0>20}.state",
        options.results_directory, options.filename_stem, options.output_count
    );
    let state = FemState::read(&path_state)?;

    // generate VTU file
    let path_vtu = format!(
        "{}/{}-{:0>20}.vtu",
        options.results_directory, options.filename_stem, options.output_count
    );
    // paraview_write_vtu(&input, &state, &path_vtu)
    Ok(())
}
