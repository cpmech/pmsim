use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use pmsim::StrError;
use structopt::StructOpt;

const DIR: &str = "/tmp/pmsim/spo_753";
const NAME: &str = "spo_753_circ_plate";

/// Command line options
#[derive(Clone, StructOpt)]
struct Options {
    /// Whether to use the arclength method or the natural method
    #[structopt(long)]
    arclength: bool,

    /// Whether to use Lagrange multipliers or static condensation for the multi-point constraints
    #[structopt(long)]
    lmm: bool,

    /// Genie for solving linear systems
    #[structopt(long, short, default_value = "mumps")]
    genie: String,

    /// Use bordering for the arclength method (only relevant if --arclength is set)
    #[structopt(long)]
    bordering: bool,

    /// Compare LMM vs SPS
    #[structopt(long)]
    lmm_vs_sps: bool,

    /// Compare bordering vs full system (only relevant if --arclength is set)
    #[structopt(long)]
    bord_vs_full: bool,
}

/// Main function
fn main() -> Result<(), StrError> {
    // first results
    let opt1 = Options::from_args();
    let mut name1 = NAME.to_string() + "_";
    name1 += &opt1.key();

    // second results
    let mut opt2 = opt1.clone();
    if opt1.lmm_vs_sps {
        opt2.lmm = !opt2.lmm;
    } else if opt1.bord_vs_full {
        opt2.bordering = !opt2.bordering;
    }
    let mut name2 = NAME.to_string() + "_";
    name2 += &opt2.key();

    // load summary and associated files
    let (post1, _) = PostProc::<2>::new(DIR, &name1)?;
    let (post2, _) = PostProc::<2>::new(DIR, &name2)?;
    let mesh = post1.mesh();
    let schema = post1.schema();

    // identify features
    let features = Features::new(&mesh, false);
    let center = features.search_point_ids(At::XY(0.0, 0.0), any_x)?[0];

    // boundaries
    let iy = schema.dof_number(center, Dof::Uy)?;

    // load results
    let nstation1 = post1.nfile();
    let nstation2 = post2.nfile();
    let mut load1 = vec![0.0; nstation1];
    let mut load2 = vec![0.0; nstation2];
    let mut deflection1 = vec![0.0; nstation1];
    let mut deflection2 = vec![0.0; nstation2];
    for index in 0..nstation1 {
        let state = post1.read_file(index)?;
        let pp = state.lambda;
        load1[index] = pp;
        deflection1[index] = -state.uu[iy];
    }
    for index in 0..nstation2 {
        let state = post2.read_file(index)?;
        let pp = state.lambda;
        load2[index] = pp;
        deflection2[index] = -state.uu[iy];
    }

    // plot
    let mut curve1 = Curve::new();
    curve1
        .set_line_style("-")
        .set_line_color("black")
        .set_marker_style(".")
        .set_marker_color("black")
        .set_label(&opt1.title())
        .draw(&deflection1, &load1);
    let mut curve2 = Curve::new();
    curve2
        .set_line_style("-")
        .set_line_color("red")
        .set_marker_style("o")
        .set_marker_size(10.0)
        .set_marker_void(true)
        .set_label(&opt2.title())
        .draw(&deflection2, &load2);
    let mut plot = Plot::new();
    let comp = if opt1.bord_vs_full {
        "bord-vs-full"
    } else if opt1.lmm_vs_sps {
        "lmm-vs-sps"
    } else {
        "no-comparison"
    };
    plot.add(&curve1)
        .add(&curve2)
        .grid_labels_legend("w (central deflection)", "P (distributed load intensity)")
        .set_figure_size_points(600.0, 250.0)
        .save(&format!("{}/{}-{}.svg", DIR, name1, comp))?;
    Ok(())
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.arclength {
            "arc".to_string()
        } else {
            "nat".to_string()
        };
        if self.lmm {
            buf += "_lmm";
        } else {
            buf += "_sps";
        }
        if self.bordering {
            buf += "_bord";
        } else {
            buf += "_full";
        }
        buf += &format!("_{}", self.genie.to_string());
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}
