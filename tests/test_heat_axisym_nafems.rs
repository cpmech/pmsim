use gemlab::prelude::*;
use plotpy::Plot;
use pmsim::{prelude::*, StrError};

fn show_mesh(mesh: &Mesh) -> Result<(), StrError> {
    let mut draw = Draw::new();
    let mut plot = Plot::new();
    draw.canvas_points.set_marker_size(3.0).set_marker_line_color("none");
    draw.cells(&mut plot, &mesh, true)?;
    draw.points(&mut plot, &mesh);
    plot.set_equal_axes(true)
        .set_figure_size_points(400.0, 600.0)
        .set_labels("x", "y")
        .save("/tmp/pmsim/mesh_heat_axisym_nafems.png")
}

#[test]
fn test_heat_axisym_nafems() -> Result<(), StrError> {
    // From Mathematica Heat Transfer Model Verification Tests
    // 2D Axisymmetric Single Equation
    // HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0001

    // geometry
    let (rin, xref, rout) = (0.02, 0.04, 0.1);
    let (ya, yb, h) = (0.04, 0.1, 0.14);

    // mesh and boundary features
    const GENERATE_MESH: bool = false;
    let mesh = if GENERATE_MESH {
        let y = &[0.0, ya, yb, h];
        let att = &[1, 1, 1];
        let (na, nb, ny) = (4, 12, &[8, 12, 8]);
        let mesh = Structured::rectangle(rin, Some(xref), rout, na, nb, y, ny, att, GeoKind::Qua9)?;
        show_mesh(&mesh)?;
        mesh.write("/tmp/pmsim/mesh_heat_axisym_nafems.dat")?;
        mesh
    } else {
        Mesh::read("data/meshes/mesh_heat_axisym_nafems.dat")?
    };

    // features
    let find = Find::new(&mesh, Some(Extract::All)); // need "All" to find reference point in the interior
    let edges_temp = find.many_edges(&[At::Y(0.0), At::Y(h), At::X(rout)], any)?;
    let edges_flux = find.edges(At::X(rin), |x| x[1] >= ya && x[1] <= yb)?;
    let ref_points = find.point_ids(At::XY(0.04, 0.04), any).unwrap();
    assert_eq!(ref_points.len(), 1);

    // reference solution
    let ref_point = ref_points[0];
    let ref_temperature = 332.97;

    // parameters, DOFs, and configuration
    let (kx, ky) = (52.0, 52.0);
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
    essential.on(&edges_temp, Ebc::T(|_| 273.15));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&edges_flux, Nbc::Qt(|_| 5e5));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // output results
    // let output = Output::new(&data);
    // output.write_vtu(&state, "/tmp/pmsim/test_heat_axisym_nafems.vtu")?;

    // check
    let eq = data.equations.eq(ref_point, Dof::T).unwrap();
    let rel_err = f64::abs(state.uu[eq] - ref_temperature) / ref_temperature;
    println!(
        "\nT = {:?}, reference = {:?}, rel_error = {:>.8} %",
        state.uu[eq],
        ref_temperature,
        rel_err * 100.0
    );
    assert!(rel_err < 1e-5);
    Ok(())
}
