use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::approx_eq;

// Arpaci's Example 3-8 on page 130 (variable conductivity)
//
// Arpaci V. S. (1966) Conduction Heat Transfer,
// Addison-Wesley, 551p
//
// TEST GOAL
//
// This tests verifies the nonlinear solver for the diffusion equation
// with a variable conductivity coefficient.
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
// Temperature T = 0 on right side @ x = L
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// Source = 5
// Variable conductivity (k = (1 + β T) kᵣ I) with kᵣ = 2
//
// NOTE
//
// The temperature at the right T = 0 (T_inf) must be zero in order to
// result in k(T_inf) = kᵣ as required by the analytical solution.

const NAME: &str = "test_heat_arpaci_nonlinear_1d";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = false;

const L: f64 = 10.0;
const SOURCE: f64 = 5.0;
const K_R: f64 = 2.0;
const BETA: f64 = 0.01;
const COEF: f64 = BETA * SOURCE * L * L / (2.0 * K_R);

// normalized analytical solution
fn normalized(x: f64) -> f64 {
    if BETA == 0.0 {
        1.0 - f64::powf(x / L, 2.0)
    } else {
        (f64::sqrt(1.0 + 2.0 * COEF * (1.0 - f64::powf(x / L, 2.0))) - 1.0) / COEF
    }
}

// analytical solution
fn analytical(x: f64) -> f64 {
    normalized(x) * SOURCE * L * L / (2.0 * K_R)
}

#[test]
fn test_heat_arpaci_nonlinear_1d() -> Result<(), StrError> {
    println!("\n################################### OLD SOLVER ###################################\n");
    run_test(false, false, false)?;
    println!("\n##################################### NATURAL ####################################\n");
    run_test(true, false, false)?; // Natural continuation
    println!("\n################################ ARCLENGTH FULL ##################################\n");
    run_test(true, true, false)?; // Pseudo-arclength continuation without bordering
    println!("\n############################# ARCLENGTH BORDERING ################################\n");
    run_test(true, true, true)?; // Pseudo-arclength continuation with bordering
    Ok(())
}

fn run_test(new_solver: bool, arclength: bool, bordering: bool) -> Result<(), StrError> {
    // mesh
    let mesh = generate_or_read_mesh(L, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let right = features.search_edges(At::X(L), any_x)?;

    // parameters
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::IsotropicLinear { kr: K_R, beta: BETA },
        source: Some(SOURCE),
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut essential = BcEssential::new();
    essential.edges(&right, Dof::Phi, 0.0); // must be zero to match analytical solution

    // natural boundary conditions
    let natural = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(false)
        .set_out_files("/tmp/pmsim", NAME, 1.0);

    // nonlinear solver configuration

    // solution
    let mut tol = 1e-13;
    let state = if new_solver {
        tol = 1e-9;
        let mut nl_config = NlConfig::new();
        nl_config
            .set_verbose(true, true, false)
            .set_h_ini(1.0)
            .set_tg_control_atol_and_rtol(0.05)
            .set_record_iterations_residuals(true);
        if arclength {
            nl_config.set_method(NlMethod::Arclength).set_bordering(bordering);
        }
        let (mut solver, mut data) = Solver::new(&mesh, &schema, &config, &essential, &natural, &mut nl_config)?;
        solver.steady(&mut data, IniDir::Pos, Stop::MaxLambda(1.0), AutoStep::Yes, None)?;
        data.get_state().clone()
    } else {
        SolverOld::solve(&mesh, &schema, &config, &essential, &natural)?
    };

    // check
    let ref_id = 0;
    let ref_x = mesh.points[ref_id].coords[0];
    let ref_eq = schema.get_eq(ref_id, Dof::Phi)?;
    let ref_tt = state.u[ref_eq];
    println!("\nT({}) = {}  ({})", ref_x, ref_tt, analytical(ref_x));
    let err = f64::abs(ref_tt - analytical(ref_x));
    println!("error = {:.5e}", err);
    approx_eq(ref_tt, analytical(ref_x), tol);

    // plot the results
    if SAVE_FIGURE {
        do_plot()
    } else {
        Ok(())
    }
}

fn do_plot() -> Result<(), StrError> {
    // get temperature values along x
    let (post, _) = PostProc::new("/tmp/pmsim", NAME)?;
    let features = Features::new(post.mesh(), false);
    let state = post.read_state(post.nstate() - 1)?;
    let (_, x_values, tt_values) = post.values_along_x(&features, &state, Dof::Phi, 0.0, any_x)?;

    // compute plot data
    let xx: Vec<_> = x_values.iter().map(|x| x / L).collect();
    let yy_num: Vec<_> = tt_values.iter().map(|tt| 2.0 * K_R * tt / (SOURCE * L * L)).collect();
    let yy_ana: Vec<_> = x_values.iter().map(|x| normalized(*x)).collect();

    // figure
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
    plot.set_title(format!("$\\beta\\;s\\;L^2\\;/\\;(2\\;k_r)$ = {:.2}", COEF).as_str())
        .grid_and_labels("$x\\;/\\;L$", "$2\\,k_r\\,T\\;/\\;(s\\,L^2)$")
        .legend()
        .save(&format!("/tmp/pmsim/{}.svg", NAME))
}

/// Generate or read mesh
fn generate_or_read_mesh(ll: f64, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[0.0, 0.0], [ll, 0.0], [ll, 1.0], [0.0, 1.0]]).unwrap();
        block.set_ndiv(&[10, 1]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(true)
            .show_cell_ids(true)
            .all(&mesh, &format!("/tmp/pmsim/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
