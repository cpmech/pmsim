use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;

// Felippa's Benchmark 14.1 (Figure 14.1) on page 14-3
//
// Felippa C, Advanced Finite Elements
//
// TEST GOAL
//
// This test verifies the axisymmetric modelling of a chick cylindrical tube
// under internal pressure. There is an analytical solution, developed for the
// plane-strain case. However, this tests employs the AXISYMMETRIC representation.
//
// MESH
//
//             Uy FIXED
//  →o------o------o------o------o
//  →|      |   .......   |      |
//  →o------o------o------o------o
//             Uy FIXED
//
// BOUNDARY CONDITIONS
//
// Fix bottom edge vertically
// Fix top edge vertically
// Distributed load Qn = -PRESSURE on left edge
//
// CONFIGURATION AND PARAMETERS
//
// Static simulation
// Young = 1000, Poisson = 0.0
// Axisymmetric
// NOTE: using 4 integration points because it gives better results with Qua8

const NAME: &str = "test_solid_felippa_thick_cylinder_axisym";

#[test]
fn test_solid_felippa_thick_cylinder_axisym() -> Result<(), StrError> {
    // Example from Felippa's A-FEM page 14-3
    const PRESSURE: f64 = 10.0;

    // mesh
    let (rin, rout, thickness) = (4.0, 10.0, 2.0);
    let mesh = generate_or_read_mesh(rin, rout, thickness, false);

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(rin), any_x)?;
    let bottom = feat.search_edges(At::Y(0.0), any_x)?;
    let top = feat.search_edges(At::Y(thickness), any_x)?;

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
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1))])?;
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
    let mut state = FemState::new(&input, &config)?;

    // run simulation
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.run(&mut state)?;

    // Felippa's Equation 14.2 on page 14-4
    let analytical_ur = |r: f64| {
        PRESSURE * (rin * rin * (1.0 + POISSON) * (rout * rout + r * r * (1.0 - 2.0 * POISSON)))
            / (YOUNG * (rout * rout - rin * rin) * r)
    };

    // check displacements
    println!("");
    let selection = feat.search_point_ids(At::Y(0.0), any_x)?;
    for p in &selection {
        let r = mesh.points[*p].coords[0];
        let eq = input.equations.eq(*p, Dof::Ux).unwrap();
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
        mesh.write(&["/tmp/pmsim/", NAME].concat()).unwrap();

        // write figure
        mesh.draw(None, &["/tmp/pmsim/", NAME, "_mesh"].concat(), |_, _| {})
            .unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&["data/meshes/", NAME, ".mesh"].concat()).unwrap()
    }
}
