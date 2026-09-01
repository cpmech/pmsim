use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::format_scientific;
use russell_sparse::Genie;
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

    /// Use bordering for the arclength method (only relevant if --arclength is set)
    #[structopt(long)]
    bordering: bool,

    #[structopt(long)]
    tol: f64,
}

/// Main function
fn main() -> Result<(), StrError> {
    // first results
    let genie1 = Genie::Cudss;
    let opt1 = Options::from_args();
    let mut name1 = NAME.to_string() + "_";
    name1 += &opt1.key(genie1);

    // second results
    let genie2 = Genie::Mumps;
    let opt2 = opt1.clone();
    let mut name2 = NAME.to_string() + "_";
    name2 += &opt2.key(genie2);

    // third results
    let genie3 = Genie::Umfpack;
    let opt3 = opt1.clone();
    let mut name3 = NAME.to_string() + "_";
    name3 += &opt3.key(genie3);

    // load summary and associated files
    let (post1, _) = PostProc::<2>::new(DIR, &name1)?;
    let (post2, _) = PostProc::<2>::new(DIR, &name2)?;
    let (post3, _) = PostProc::<2>::new(DIR, &name3)?;
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
    let nstation3 = post3.nfile();
    let mut load1 = vec![0.0; nstation1];
    let mut load2 = vec![0.0; nstation2];
    let mut load3 = vec![0.0; nstation3];
    let mut deflection1 = vec![0.0; nstation1];
    let mut deflection2 = vec![0.0; nstation2];
    let mut deflection3 = vec![0.0; nstation3];
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
    for index in 0..nstation3 {
        let state = post3.read_file(index)?;
        let pp = state.lambda;
        load3[index] = pp;
        deflection3[index] = -state.uu[iy];
    }

    // plot
    let mut curve1 = Curve::new();
    curve1
        .set_line_style("-")
        .set_line_color("red")
        .set_marker_style("o")
        .set_marker_color("red")
        .set_label(&genie1.to_string().to_uppercase())
        .draw(&deflection1, &load1);
    let mut curve2 = Curve::new();
    curve2
        .set_line_style("-")
        .set_line_color("green")
        .set_marker_style("^")
        .set_marker_size(9.0)
        .set_marker_void(true)
        .set_label(&genie2.to_string().to_uppercase())
        .draw(&deflection2, &load2);
    let mut curve3 = Curve::new();
    curve3
        .set_line_style("-")
        .set_line_color("blue")
        .set_marker_style("s")
        .set_marker_color("blue")
        .set_marker_size(7.0)
        .set_marker_void(true)
        .set_label(&genie3.to_string().to_uppercase())
        .draw(&deflection3, &load3);
    let mut plot = Plot::new();
    let str_tol = format_scientific(opt1.tol, 7, 1);
    let title = format!("TOL RDIFF MIN = {}", str_tol);
    plot.add(&curve1)
        .add(&curve2)
        .add(&curve3)
        .set_title(&title)
        .grid_labels_legend("w (central deflection)", "P (distributed load intensity)")
        .save(&format!("{}/{}-mumps-umfpack-{}.svg", DIR, name1, str_tol))?;
    Ok(())
}

impl Options {
    fn key(&self, genie: Genie) -> String {
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
        buf += &format!("_{}", genie.to_string());
        buf
    }
}
