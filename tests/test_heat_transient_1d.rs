use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};
use russell_lab::math::{erfc, PI};

#[test]
fn test_heat_transient_1d() -> Result<(), StrError> {
    // mesh
    const GENERATE_MESH: bool = false;
    let mesh = if GENERATE_MESH {
        let mut block = Block::new(&[[0.0, 0.0], [20.0, 0.0], [20.0, 1.0], [0.0, 1.0]])?;
        block.set_ndiv(&[10, 1])?;
        let mesh = block.subdivide(GeoKind::Qua8)?;
        mesh.write("/tmp/pmsim/mesh_heat_transient_1d_qua8.dat")?;
        draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_transient_1d_qua8.svg")?;
        mesh
    } else {
        // let mesh = Mesh::read("data/meshes/mesh_heat_transient_1d.dat")?;
        Mesh::read("data/meshes/mesh_heat_transient_1d_qua8.dat")?
    };

    // features
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(0.0), any)?;

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
        find.point_ids(At::X(0.0), any).unwrap(),
        find.point_ids(At::X(1.0), any).unwrap(),
        find.point_ids(At::X(2.0), any).unwrap(),
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
    Ok(())
}
