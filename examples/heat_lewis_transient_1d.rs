use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::{prelude::*, StrError};
use russell_lab::math::{erfc, PI};
use russell_lab::prelude::*;

const FILENAME_KEY: &'static str = "test_heat_lewis_transient_1d";

// Lewis' Example 6.4.2 on page 159
//
// Lewis R, Nithiarasu P, and Seetharamu KN (2004) Fundamentals of the
// Finite Element Method for Heat and Fluid Flow, Wiley, 341p
//
// TEST GOAL
//
// This test verifies the transient diffusion in 1D
//
// MESH
//
// o-----------------------------------------------------------o
// |    |    |    |    |    |    |    |    |    |    .....     | h = 1
// o-----------------------------------------------------------o
//                      <-  L = 20 ->
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Flux Qt = 1 on left side @ x = 0
//
// CONFIGURATION AND PARAMETERS
//
// Transient simulation
// No source
// Constant conductivity kx = ky = 1
// Coefficient ρ = 1

fn main() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&FilePath::mesh(FILENAME_KEY, false))?;

    // features
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(0.0), any_x)?;

    // parameters, DOFs, and configuration
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::Constant {
            kx: 1.0,
            ky: 1.0,
            kz: 1.0,
        },
        source: None,
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let mut config = Config::new();
    let t_fin = 1.0;
    config.transient = true;
    config.control.t_fin = t_fin;

    // essential boundary conditions
    let essential = Essential::new();

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left, Nbc::Qt(|_| 1.0));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // check
    let analytical = |t: f64, x: f64| {
        2.0 * f64::sqrt(t / PI)
            * (f64::exp(-x * x / (4.0 * t)) - (x / 2.0) * f64::sqrt(PI / t) * erfc(x / (2.0 * f64::sqrt(t))))
    };
    let selected = vec![
        find.point_ids(At::X(0.0), any_x).unwrap(),
        find.point_ids(At::X(1.0), any_x).unwrap(),
        find.point_ids(At::X(2.0), any_x).unwrap(),
    ]
    .concat();
    println!("");
    for p in &selected {
        let x = data.mesh.points[*p].coords[0];
        let eq = data.equations.eq(*p, Dof::T).unwrap();
        let tt = state.uu[eq];
        let diff = f64::abs(tt - analytical(state.t, x));
        println!("point = {}, x = {:.2}, T = {:.6}, diff = {:.4e}", p, x, tt, diff);
        assert!(diff < 3e-2);
    }

    // plot results
    if false {
        // compute analytical solution
        let xx_ana = Vector::linspace(0.0, 2.0, 11)?;
        let tt_ana = xx_ana.get_mapped(|x| analytical(t_fin, x));

        // get temperature values along x
        let post = PostProc::new(&mesh, &find, &data, &state);
        let (_, xx_num, tt_num) = post.values_along_x(Dof::T, 0.0, |x| x[0] <= 2.0)?;

        // plot
        let mut curve_ana = Curve::new();
        let mut curve_num = Curve::new();
        curve_ana.draw(&xx_ana, &tt_ana);
        curve_num
            .set_line_color("#cd0000")
            .set_line_style("None")
            .set_marker_style("+");
        curve_num.draw(&xx_num, &tt_num);
        let mut plot = Plot::new();
        plot.add(&curve_ana).add(&curve_num);
        plot.grid_and_labels("x", "T")
            .set_yrange(0.0, 1.2)
            .legend()
            .save("/tmp/pmsim/test_heat_transient_1d.svg")?;
    }
    Ok(())
}