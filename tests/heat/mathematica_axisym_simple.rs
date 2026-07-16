use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};

// From Mathematica Heat Transfer Model Verification Tests
// (HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0001)
//
// 2D Axisymmetric Single Equation
//
// https://reference.wolfram.com/language/PDEModels/tutorial/HeatTransfer/HeatTransferVerificationTests.html
//
// TEST GOAL
//
// This test verifies the steady heat equation in 1D with prescribed flux
//
// MESH
//
//   →→ ---------------------
//   →→ |    |    |    |    |  h
//   →→ ---------------------
//     1.0                 2.0
//     rin                rout
//
// INITIAL CONDITIONS
//
// Temperature T = 0 at all points
//
// BOUNDARY CONDITIONS
//
// Temperature T = 10.0 on the right edge
// Inward Flux Qt = -100.0 on the left edge
//
// CONFIGURATION AND PARAMETERS
//
// Steady simulation
// No source
// Constant conductivity kx = ky = 10.0

const NAME: &str = "mathematica_axisym_simple";
const GENERATE_MESH: bool = false;

#[test]
fn mathematica_axisym_simple() -> Result<(), StrError> {
    // geometry
    let (rin, rout, h) = (1.0, 2.0, 0.1);

    // mesh
    let mesh = generate_or_read_mesh(rin, rout, h, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(rin), any_x)?;
    let right = features.search_edges(At::X(rout), any_x)?;

    // parameters
    let (kx, ky) = (10.0, 10.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&right, Dof::Phi, 10.0);

    // natural boundary conditions
    let mut nbc = BcNatural::new();
    nbc.edges(&left, Nbc::Qt, -100.0); // inward flux

    // configuration
    let mut config = Config::new(&mesh);
    config.axisymmetric();

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;
    let state = data.state();

    // check
    let analytical = |r: f64| 10.0 * (1.0 - f64::ln(r / 2.0));
    for point in &mesh.points {
        let x = point.coords[0];
        let i = schema.dof_number(point.id, Dof::Phi)?;
        let tt = state.uu[i];
        let diff = f64::abs(tt - analytical(x));
        // println!("point = {}, x = {:.2}, T = {:.6}, diff = {:.4e}", point.id, x, tt, diff);
        assert!(diff < 1e-5);
    }
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(rin: f64, rout: f64, h: f64, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[rin, 0.0], [rout, 0.0], [rout, h], [rin, h]]).unwrap();
        block.set_ndiv(&[10, 1]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua9).unwrap();

        // draw figure
        let mut draw = Draw::new();
        draw.show_point_ids(true)
            .show_cell_ids(true)
            .set_range_2d(0.95, 2.05, -0.05, 0.15)
            .set_size(600.0, 100.0)
            .all(&mesh, &format!("/tmp/pmsim/heat/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/heat/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
