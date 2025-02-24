use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use russell_lab::*;

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
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&right, Dof::Phi, 0.0); // must be zero to match analytical solution

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let config = Config::new(&mesh);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check
    let ref_id = 0;
    let ref_x = mesh.points[ref_id].coords[0];
    let ref_eq = base.dofs.eq(ref_id, Dof::Phi)?;
    let ref_tt = state.u[ref_eq];
    println!("\nT({}) = {}  ({})", ref_x, ref_tt, analytical(ref_x));
    approx_eq(ref_tt, analytical(ref_x), 1e-13);

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
    let state = post.read_state(post.n_state() - 1)?;
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
        let mut fig = Figure::new();
        fig.show_point_ids(true)
            .show_cell_ids(true)
            .draw(&mesh, &format!("/tmp/pmsim/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
