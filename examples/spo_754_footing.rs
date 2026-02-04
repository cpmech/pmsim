use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{Stopwatch, Vector};

const NAME: &str = "spo_754_footing";

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

// loading factors
const LAMBDAS: [f64; 16] = [
    0.0,   //  0
    0.01,  //  1
    0.015, //  2
    0.02,  //  3
    0.025, //  4
    0.035, //  5
    0.045, //  6
    0.055, //  7
    0.065, //  8
    0.075, //  9
    0.08,  // 10
    0.085, // 10b (something happens that needs this extra increment)
    0.09,  // 11
    0.11,  // 12
    0.14,  // 13
    0.2,   // 14
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
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&footing, Dof::Uy, -1.0);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config.out_files("/tmp/pmsim", NAME);

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(true, true, true)
        .set_record_iterations_residuals(true);

    // simulator
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // run simulation
    let stop = Stop::Steps(LAMBDAS.len() - 1);
    let list = Vector::from(&LAMBDAS).get_differences();
    let dll = DeltaLambda::list(list.as_data());
    sim.steady(&mut data, IniDir::Pos, stop, dll)?;

    // stop stopwatch
    sw.stop();
    println!("\nelapsed time = {}\n", sw);
    Ok(())
}
