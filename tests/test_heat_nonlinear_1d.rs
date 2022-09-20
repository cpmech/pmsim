use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::{prelude::*, StrError};
use russell_chk::approx_eq;

//
//
//
// MESH
//
// o-----------------------------------------------------------o
// |    |    |    |    |    |    |    |    |    |    .....     | h = 1
// o-----------------------------------------------------------o
//                      <-  L = 10 ->
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Temperature T on right side @ x = L
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// Source = 5
// Variable conductivity (k = (1 + β T) kᵣ I) with kᵣ = 2

#[test]
fn test_heat_nonlinear_1d() -> Result<(), StrError> {
    // constants
    const L: f64 = 10.0;
    const SOURCE: f64 = 5.0;
    const K_R: f64 = 2.0;
    const BETA: f64 = 0.01;

    // mesh
    const GENERATE_MESH: bool = false;
    let mesh = if GENERATE_MESH {
        let mut block = Block::new(&[[0.0, 0.0], [L, 0.0], [L, 1.0], [0.0, 1.0]])?;
        block.set_ndiv(&[10, 1])?;
        let mesh = block.subdivide(GeoKind::Qua4)?;
        mesh.write("/tmp/pmsim/mesh_heat_nonlinear_1d.dat")?;
        draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_nonlinear_1d.svg")?;
        mesh
    } else {
        Mesh::read("data/meshes/mesh_heat_nonlinear_1d.dat")?
    };

    // features
    let find = Find::new(&mesh, None);
    let right = find.edges(At::X(L), any_x)?;

    // parameters, DOFs, and configuration
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::IsotropicLinear { kr: K_R, beta: BETA },
        source: Some(SOURCE),
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&right, Ebc::T(|_| 0.0)); // must be zero to match analytical solution

    // natural boundary conditions
    let natural = Natural::new();

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // analytical solution
    let coef = BETA * SOURCE * L * L / (2.0 * K_R);
    let normalized = |x: f64| {
        if BETA == 0.0 {
            1.0 - f64::powf(x / L, 2.0)
        } else {
            (f64::sqrt(1.0 + 2.0 * coef * (1.0 - f64::powf(x / L, 2.0))) - 1.0) / coef
        }
    };
    let analytical = |x: f64| normalized(x) * SOURCE * L * L / (2.0 * K_R);

    // check
    let ref_id = 0;
    let ref_x = mesh.points[ref_id].coords[0];
    let ref_eq = data.equations.eq(ref_id, Dof::T)?;
    let ref_tt = state.uu[ref_eq];
    println!("\nT({}) = {}  ({})", ref_x, ref_tt, analytical(ref_x));
    approx_eq(ref_tt, analytical(ref_x), 1e-13);

    // plot results
    if false {
        // get temperature values along x
        let post = PostProc::new(&mesh, &find, &data, &state);
        let (_, x_values, tt_values) = post.values_along_x(Dof::T, 0.0, any_x)?;

        // compute plot data
        let xx: Vec<_> = x_values.iter().map(|x| x / L).collect();
        let yy_num: Vec<_> = tt_values.iter().map(|tt| 2.0 * K_R * tt / (SOURCE * L * L)).collect();
        let yy_ana: Vec<_> = x_values.iter().map(|x| normalized(*x)).collect();

        // plot
        let mut curve_num = Curve::new();
        let mut curve_ana = Curve::new();
        curve_num
            .set_line_color("#cd0000")
            .set_line_style("None")
            .set_marker_style("+");
        curve_num.draw(&xx, &yy_num);
        curve_ana.draw(&xx, &yy_ana);
        let mut plot = Plot::new();
        plot.add(&curve_ana);
        plot.add(&curve_num);
        plot.set_title(format!("$\\beta\\;s\\;L^2\\;/\\;(2\\;k_r)$ = {:.2}", coef).as_str())
            .grid_and_labels("$x\\;/\\;L$", "$2\\,k_r\\,T\\;/\\;(s\\,L^2)$")
            .legend()
            .save("/tmp/pmsim/test_heat_nonlinear_1d.svg")?;
    }
    Ok(())
}
