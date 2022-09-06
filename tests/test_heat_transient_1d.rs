use gemlab::prelude::*;
use pmsim::{fem::sim_transient, prelude::*, StrError};
use russell_lab::math::{erfc, PI};

#[test]
fn test_heat_transient_1d() -> Result<(), StrError> {
    // mesh and boundary features
    // let mesh = Mesh::read("data/meshes/mesh_heat_transient_1d.dat")?;
    let mesh = Mesh::read("data/meshes/mesh_heat_transient_1d_qua8.dat")?;
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(0.0))?;

    // parameters, DOFs, and configuration
    let p1 = ParamDiffusion {
        rho: 1.0,
        kx: 1.0,
        ky: 1.0,
        kz: 1.0,
        source: None,
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let mut config = Config::new();
    config.transient = true;

    // essential boundary conditions
    let essential = Essential::new();

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left, Nbc::Qt(|_| 1.0));

    // interior elements
    let mut interior_elements = InteriorElementVec::new(&data, &config)?;

    // boundary elements
    let mut boundary_elements = BoundaryElementVec::new(&data, &config, &natural)?;

    // simulation state
    let mut state = State::new(&data, &config, &essential)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).unwrap();

    // run simulation
    sim_transient(
        &mut interior_elements,
        &mut boundary_elements,
        &mut state,
        &mut lin_sys,
        &config,
    )?;

    // check
    let analytical = |t: f64, x: f64| {
        2.0 * f64::sqrt(t / PI)
            * (f64::exp(-x * x / (4.0 * t)) - (x / 2.0) * f64::sqrt(PI / t) * erfc(x / (2.0 * f64::sqrt(t))))
    };
    let selected = vec![
        find.point_ids(At::X(0.0)).unwrap(),
        find.point_ids(At::X(1.0)).unwrap(),
        find.point_ids(At::X(2.0)).unwrap(),
    ]
    .concat();
    println!("");
    for p in &selected {
        let x = data.mesh.points[*p].coords[0];
        let eq = data.equations.eq(*p, Dof::T).unwrap();
        let tt = state.uu[eq];
        let diff = f64::abs(tt - analytical(state.t, x));
        println!("point = {}, x = {:.2}, tt = {:.6}, diff = {:.4e}", p, x, tt, diff);
        assert!(diff < 3e-2);
    }
    Ok(())
}
