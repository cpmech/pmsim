use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};

#[test]
fn test_heat_axisym_simple() -> Result<(), StrError> {
    // From Mathematica Heat Transfer Model Verification Tests
    // 2D Axisymmetric Single Equation
    // HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0001

    // mesh and boundary features
    let (rin, rout, h) = (1.0, 2.0, 0.1);
    let mut block = Block::new(&[[rin, 0.0], [rout, 0.0], [rout, h], [rin, h]])?;
    block.set_ndiv(&[10, 1])?;
    let mesh = block.subdivide(GeoKind::Qua17)?;
    // draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_axisym_simple_2d.svg")?;
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(rin), any)?;
    let right = find.edges(At::X(rout), any)?;

    // parameters, DOFs, and configuration
    let (kx, ky) = (10.0, 10.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        kx,
        ky,
        kz: 0.0,
        source: None,
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let mut config = Config::new();
    config.axisymmetric = true;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&right, Ebc::T(|_| 10.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left, Nbc::Qt(|_| 100.0));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;
    println!("{}", state.uu);

    // check
    let analytical = |r: f64| 10.0 * (1.0 - f64::ln(r / 2.0));
    for point in &mesh.points {
        let x = point.coords[0];
        let eq = data.equations.eq(point.id, Dof::T).unwrap();
        let tt = state.uu[eq];
        let diff = f64::abs(tt - analytical(x));
        println!("point = {}, x = {:.2}, T = {:.6}, diff = {:.4e}", point.id, x, tt, diff);
        assert!(diff < 1e-7);
    }
    Ok(())
}
