use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::Stopwatch;

const NAME: &str = "spo_754_footing";

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

// displacement control
const UY: [f64; 15] = [
    0.0,    //
    -0.01,  //
    -0.015, //
    -0.02,  //
    -0.025, //
    -0.035, //
    -0.045, //
    -0.055, //
    -0.065, //
    -0.075, //
    -0.08,  //
    -0.09,  //
    -0.11,  //
    -0.14,  //
    -0.2,   //
];

pub fn main() -> Result<(), StrError> {
    // start stopwatch
    let mut sw = Stopwatch::new();

    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(500.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let footing = features.search_edges(At::Y(500.0), |x| x[0] <= 50.0)?;

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            // stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&footing, Dof::Uy, |t| UY[t as usize]);

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_lagrange_mult_method(true)
        .set_steady(UY.len() - 1)
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_max_iterations(20);

    // solution
    SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?;

    // stop stopwatch
    sw.stop();
    println!("\nelapsed time = {}\n", sw);
    Ok(())
}
