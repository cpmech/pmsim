use gemlab::mesh::Mesh;
use pmsim::fem::{FemState, FileIo};
use pmsim::StrError;
use std::fmt::Write;
use std::fs::File;
use std::io::Write as IoWrite;
use structopt::StructOpt;

/// Command line options
#[derive(StructOpt, Debug)]
#[structopt(
    name = "pmsim_to_paraview",
    about = "Generates VTU files for visualization with Paraview"
)]
struct Options {
    /// Directory containing the result files
    results_dir: String,

    /// Filename stem
    filename_stem: String,

    /// Index of the output ("-1" means all)
    #[structopt(short = "i", long, default_value = "-1")]
    output_index: i32,
}

fn main() -> Result<(), StrError> {
    // parse options
    let options = Options::from_args();

    // load FileIo
    let file_io = FileIo::read_json(&format!(
        "{}/{}-summary.json",
        options.results_dir, options.filename_stem
    ))?;

    /*
    // selected indices
    let indices = if options.output_index < 0 {
        file_io.indices
    } else {
        vec![options.output_index as usize]
    };

    // load mesh
    let mesh = Mesh::read_json(&format!("{}/{}-mesh.json", options.results_dir, options.filename_stem))?;

    // initiate PVD fle
    let mut pvd = String::new();
    write!(&mut pvd, "<?xml version=\"1.0\"?>\n<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\">\n<Collection>\n").unwrap();

    // generate VTU file
    for index in &indices {
        let state = FemState::read_json(&format!(
            "{}/{}-{:0>20}.json",
            options.results_dir, options.filename_stem, index
        ))?;
        let vtu_fn = file_io.path_vtu(*index);
        write!(
            &mut pvd,
            "<DataSet timestep=\"{:?}\" file=\"{}\" />\n",
            file_io.times[*index], vtu_fn
        )
        .unwrap();
        file_io.write_vtu(&mesh, &state, *index)?;
    }

    // generate PVD file
    write!(&mut pvd, "</Collection>\n</VTKFile>").unwrap();
    let pvd_fn = format!("{}/{}.pvd", options.results_dir, options.filename_stem);
    let mut file = File::create(&pvd_fn).map_err(|_| "cannot create file")?;
    file.write_all(pvd.as_bytes()).map_err(|_| "cannot write file")?;
    file.sync_all().map_err(|_| "cannot sync file")?;

    // message
    let thin_line = format!("{:─^1$}", "", pvd_fn.len());
    println!("\n\n{}", thin_line);
    println!("VTU files generated; the PVD file is:");
    println!("{}", pvd_fn);
    println!("{}\n\n", thin_line);
    */
    Ok(())
}
