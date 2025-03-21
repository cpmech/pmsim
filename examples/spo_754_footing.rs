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
const UY: [f64; 14] = [
    -0.01,  // stage =  0
    -0.015, // stage =  1
    -0.02,  // stage =  2
    -0.025, // stage =  3
    -0.035, // stage =  4
    -0.045, // stage =  5
    -0.055, // stage =  6
    -0.065, // stage =  7
    -0.075, // stage =  8
    -0.08,  // stage =  9
    -0.09,  // stage = 10
    -0.11,  // stage = 11
    -0.14,  // stage = 12
    -0.2,   // stage = 13
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
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&footing, Dof::Uy, 1.0, |stage, _| UY[stage]);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_out_files("/tmp/pmsim", NAME, 1.0)
        .set_lagrange_mult_method(true)
        .set_steady(UY.len())
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &base, &config)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // stop stopwatch
    sw.stop();
    println!("\nelapsed time = {}\n", sw);
    Ok(())
}
