use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::prelude::*;
use russell_lab::math::{erfc, PI};
use russell_lab::*;

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

const NAME: &str = "test_heat_lewis_transient_1d";

const SAVE_FIGURE: bool = false;

#[test]
fn test_heat_lewis_transient_1d() -> Result<(), StrError> {
    // mesh
    let mesh = generate_or_read_mesh(false);

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;

    // parameters
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant {
            kx: 1.0,
            ky: 1.0,
            kz: 1.0,
        },
        source: None,
        ngauss: None,
    };
    let fem = FemMesh::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let essential = Essential::new();

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&left, Nbc::Qt, 1.0);

    // configuration
    let mut config = Config::new(&mesh);
    let t_fin = 1.0;
    config.set_transient(true).set_dt(|_| 0.1).set_t_fin(t_fin);

    // FEM state
    let mut state = FemState::new(&fem, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&fem, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // check
    let analytical = |t: f64, x: f64| {
        2.0 * f64::sqrt(t / PI)
            * (f64::exp(-x * x / (4.0 * t)) - (x / 2.0) * f64::sqrt(PI / t) * erfc(x / (2.0 * f64::sqrt(t))))
    };
    let selected = vec![
        features.search_point_ids(At::X(0.0), any_x).unwrap(),
        features.search_point_ids(At::X(1.0), any_x).unwrap(),
        features.search_point_ids(At::X(2.0), any_x).unwrap(),
    ]
    .concat();
    println!("");
    for p in &selected {
        let x = fem.mesh.points[*p].coords[0];
        let eq = fem.equations.eq(*p, Dof::T).unwrap();
        let tt = state.uu[eq];
        let diff = f64::abs(tt - analytical(state.t, x));
        println!("point = {}, x = {:.2}, T = {:.6}, diff = {:.4e}", p, x, tt, diff);
        assert!(diff < 3e-2);
    }

    // plot
    if SAVE_FIGURE {
        // compute analytical solution
        let xx_ana = Vector::linspace(0.0, 2.0, 11)?;
        let tt_ana = xx_ana.get_mapped(|x| analytical(t_fin, x));

        // get temperature values along x
        let post = PostProc::new(&fem, &config);
        let (_, xx_num, tt_num) = post.values_along_x(&features, &state, Dof::T, 0.0, |x| x[0] <= 2.0)?;

        // plot
        let mut curve_ana = Curve::new();
        let mut curve_num = Curve::new();
        curve_ana.draw(xx_ana.as_data(), tt_ana.as_data());
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
            .save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[0.0, 0.0], [20.0, 0.0], [20.0, 1.0], [0.0, 1.0]]).unwrap();
        block.set_ndiv(&[10, 1]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();

        // write mesh
        mesh.write_json(&format!("{}/{}.json", DEFAULT_TEST_DIR, NAME)).unwrap();

        // write figure
        mesh.draw(None, &format!("{}/{}_mesh.svg", DEFAULT_TEST_DIR, NAME), |_, _| {})
            .unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read_json(&format!("data/meshes/{}.json", NAME)).unwrap()
    }
}
