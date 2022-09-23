use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;

const FILENAME_KEY: &'static str = "test_solid_felippa_thick_cylinder_axisym";

#[test]
fn test_solid_felippa_thick_cylinder_axisym() -> Result<(), StrError> {
    // Example from Felippa's A-FEM page 13-5
    const PRESSURE: f64 = 10.0;

    // mesh
    let (rin, rout, thickness) = (4.0, 10.0, 2.0);
    let mesh = generate_or_read_mesh(rin, rout, thickness, false);

    // features
    let find = Find::new(&mesh, None);
    let left = find.edges(At::X(rin), any_x)?;
    let bottom = find.edges(At::Y(0.0), any_x)?;
    let top = find.edges(At::Y(thickness), any_x)?;

    const YOUNG: f64 = 1000.0;
    const POISSON: f64 = 0.0;
    // parameters, DOFs, and configuration
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(p1))])?;
    let mut config = Config::new();
    config.axisymmetric = true;
    config.n_integ_point.insert(1, 4); // reduced integration => better results

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&bottom, Ebc::Uy(|_| 0.0)).on(&top, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left, Nbc::Qn(|_| -PRESSURE));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // Felippa's Equation 14.2 on page 14-4
    let analytical_ur = |r: f64| {
        PRESSURE * (rin * rin * (1.0 + POISSON) * (rout * rout + r * r * (1.0 - 2.0 * POISSON)))
            / (YOUNG * (rout * rout - rin * rin) * r)
    };

    // check displacements
    println!("");
    let selection = find.point_ids(At::Y(0.0), any_x)?;
    for p in &selection {
        let r = mesh.points[*p].coords[0];
        let eq = data.equations.eq(*p, Dof::Ux).unwrap();
        let ux = state.uu[eq];
        let diff = f64::abs(ux - analytical_ur(r));
        println!("point = {}, r = {:?}, Ux = {:?}, diff = {:?}", p, r, ux, diff);
        assert!(diff < 1e-15);
    }
    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(rin: f64, rout: f64, thickness: f64, generate: bool) -> Mesh {
    if generate {
        // generate mesh
        let mut block = Block::new(&[[rin, 0.0], [rout, 0.0], [rout, thickness], [rin, thickness]]).unwrap();
        block.set_ndiv(&[2, 1]).unwrap();
        let mesh = block.subdivide(GeoKind::Qua8).unwrap();

        // write mesh
        mesh.write(&FilePath::mesh(FILENAME_KEY, true)).unwrap();

        // write figure
        draw_mesh(&mesh, true, &FilePath::svg_suffix(FILENAME_KEY, "_mesh", true)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&FilePath::mesh(FILENAME_KEY, false)).unwrap()
    }
}