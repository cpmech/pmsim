use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;

// NOTE: the results from Kythe's book are WRONG
#[test]
fn test_integ_heat_axisym_kythe() -> Result<(), StrError> {
    // mesh
    #[rustfmt::skip]
    let mesh = Mesh {
        ndim: 2,
        points: vec![
            Point { id: 0, coords: vec![0.0,  0.0 ] },
            Point { id: 1, coords: vec![0.04, 0.0 ] },
            Point { id: 2, coords: vec![0.04, 0.08] },
            Point { id: 3, coords: vec![0.0,  0.08] },
        ],
        cells: vec![
            Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 2, 3] },
            Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 2] },
        ],
    };
    // draw_mesh(&mesh, true, "/tmp/pmsim/mesh_integ_heat_axisym_kythe.svg")?;
    let find = Find::new(&mesh, None);
    let bottom = find.edges(At::Y(0.0))?;
    let top = find.edges(At::Y(0.08))?;
    let right = find.edges(At::X(0.04))?;

    // parameters, DOFs, and configuration
    let (kx, ky) = (0.256, 0.256);
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
    essential.on(&bottom, Ebc::T(|_| 121.111));

    // natural boundary conditions
    let convection = Nbc::Cv(35.6, |_| 121.111);
    let mut natural = Natural::new();
    natural.on(&top, convection).on(&right, convection);

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;
    println!("{}", state.uu);

    Ok(())
}
