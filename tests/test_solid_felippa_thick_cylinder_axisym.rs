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
const GENERATE_MESH: bool = false;

#[test]
fn test_solid_felippa_thick_cylinder_axisym() -> Result<(), StrError> {
    // Example from Felippa's A-FEM page 14-3
    const PRESSURE: f64 = 10.0;

    // mesh
    let (rin, rout, thickness) = (4.0, 10.0, 2.0);
    let mesh = generate_or_read_mesh(rin, rout, thickness, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(rin), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(thickness), any_x)?;

    const YOUNG: f64 = 1000.0;
    const POISSON: f64 = 0.0;
    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
        ngauss: Some(4), // reduced integration => better results
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&bottom, Dof::Uy, 0.0).edges(&top, Dof::Uy, 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges(&left, Nbc::Qn, -PRESSURE);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_axisymmetric();

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // Felippa's Equation 14.2 on page 14-4
    let analytical_ur = |r: f64| {
        PRESSURE * (rin * rin * (1.0 + POISSON) * (rout * rout + r * r * (1.0 - 2.0 * POISSON)))
            / (YOUNG * (rout * rout - rin * rin) * r)
    };

    // check displacements
    println!("");
    let selection = features.search_point_ids(At::Y(0.0), any_x)?;
    for p in &selection {
        let r = mesh.points[*p].coords[0];
        let eq = base.dofs.eq(*p, Dof::Ux).unwrap();
        let ux = state.u[eq];
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

        // draw figure
        let mut fig = Figure::new();
        fig.show_point_ids(true)
            .show_cell_ids(true)
            .draw(&mesh, &format!("/tmp/pmsim/mesh_{}.svg", NAME))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}.msh", NAME)).unwrap();
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}.msh", NAME)).unwrap()
    }
}
