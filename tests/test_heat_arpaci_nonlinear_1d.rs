use gemlab::prelude::*;
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

#[test]
fn test_heat_arpaci_nonlinear_1d() -> Result<(), StrError> {
    // constants
    const L: f64 = 10.0;
    const SOURCE: f64 = 5.0;
    const K_R: f64 = 2.0;
    const BETA: f64 = 0.01;

    // mesh
    let mesh = generate_or_read_mesh(L, false);

    // features
    let feat = Features::new(&mesh, false);
    let right = feat.search_edges(At::X(L), any_x)?;

    // input data
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::IsotropicLinear { kr: K_R, beta: BETA },
        source: Some(SOURCE),
    };
    let input = FemInput::new(&mesh, [(1, Element::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&right, Ebc::T(|_| 0.0)); // must be zero to match analytical solution

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let config = Config::new();

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

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
    let ref_eq = input.equations.eq(ref_id, Dof::T)?;
    let ref_tt = state.uu[ref_eq];
    println!("\nT({}) = {}  ({})", ref_x, ref_tt, analytical(ref_x));
    approx_eq(ref_tt, analytical(ref_x), 1e-13);
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(ll: f64, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[0.0, 0.0], [ll, 0.0], [ll, 1.0], [0.0, 1.0]]).unwrap();
        block.set_ndiv(&[10, 1]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua4).unwrap();

        // write mesh
        mesh.write_json(&format!("{}/{}.json", DEFAULT_TEST_DIR, NAME)).unwrap();

        // write figure
        mesh.draw(None, &format!("{}/{}.svg", DEFAULT_TEST_DIR, NAME), |_, _| {})
            .unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read_json(&format!("data/meshes/{}.json", NAME)).unwrap()
    }
}
