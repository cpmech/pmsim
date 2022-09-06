use gemlab::prelude::*;
use pmsim::{fem::sim_transient, prelude::*, StrError};

#[test]
fn test_heat_transient_1d() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = Mesh::read("data/meshes/mesh_heat_transient_1d.dat")?;
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
    config.control.n_max_time_steps = 2;

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
    )
}
