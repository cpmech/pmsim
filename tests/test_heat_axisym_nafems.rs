#![allow(unused)]

use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};
use russell_chk::approx_eq;

#[test]
fn test_heat_axisym_nafems() -> Result<(), StrError> {
    // From Mathematica Heat Transfer Model Verification Tests
    // 2D Axisymmetric Single Equation
    // HeatTransfer-FEM-Stationary-2DAxisym-Single-HeatTransfer-0001

    // geometry
    let (rin, rout, h, ya, yb) = (0.02, 0.1, 0.14, 0.04, 0.1);

    // mesh and boundary features
    const GENERATE_MESH: bool = true;
    const FINE_MESH: bool = true;
    let mesh = if GENERATE_MESH {
        let mut block1 = Block::new(&[[rin, 0.0], [rout, 0.0], [rout, ya], [rin, ya]])?;
        let mut block2 = Block::new(&[[rin, ya], [rout, ya], [rout, yb], [rin, yb]])?;
        let mut block3 = Block::new(&[[rin, yb], [rout, yb], [rout, h], [rin, h]])?;
        if FINE_MESH {
            // block1.set_ndiv(&[32, 16])?;
            // block2.set_ndiv(&[32, 24])?;
            // block3.set_ndiv(&[32, 16])?;
            block1.set_ndiv(&[16, 8])?;
            block2.set_ndiv(&[16, 12])?;
            block3.set_ndiv(&[16, 8])?;
            // block1.set_ndiv(&[8, 4])?;
            // block2.set_ndiv(&[8, 6])?;
            // block3.set_ndiv(&[8, 4])?;
        } else {
            block1.set_ndiv(&[4, 2])?;
            block2.set_ndiv(&[4, 3])?;
            block3.set_ndiv(&[4, 2])?;
        }
        let mesh1 = block1.subdivide(GeoKind::Qua9)?;
        let mesh2 = block2.subdivide(GeoKind::Qua9)?;
        let mesh3 = block3.subdivide(GeoKind::Qua9)?;
        let mesh = join_meshes(&[&mesh1, &mesh2, &mesh3])?;
        // draw_mesh(&mesh, false, "/tmp/pmsim/mesh_heat_axisym_nafems.svg")?;
        // draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_axisym_nafems.svg")?;
        // mesh.write("/tmp/pmsim/mesh_heat_axisym_nafems.dat")?;
        mesh
    } else {
        Mesh::read("data/meshes/mesh_heat_axisym_nafems.dat")?
    };

    // features
    let find = Find::new(&mesh, Some(Extract::All)); // need "All" to find reference point
    let bot = find.edges(At::Y(0.0), any)?;
    let top = find.edges(At::Y(h), any)?;
    let left = find.edges(At::X(rin), any)?;
    let right = find.edges(At::X(rout), any)?;
    let left_flux: Vec<_> = left
        .iter()
        .filter(|&&feature| {
            feature.points.iter().fold(true, |belong, p| {
                let y = mesh.points[*p].coords[1];
                if y < ya || y > yb {
                    false
                } else {
                    belong
                }
            })
        })
        .copied()
        .collect();
    let ref_points = find.point_ids(At::XY(0.04, 0.04), any).unwrap();
    println!("ref_points: {:?}", ref_points);
    println!("left: {:?}", left.iter().map(|f| &f.points).collect::<Vec<_>>());
    println!(
        "left_flux: {:?}\n",
        left_flux.iter().map(|f| &f.points).collect::<Vec<_>>()
    );
    assert_eq!(ref_points.len(), 1);
    println!("bot: {:?}", bot.iter().map(|f| &f.points).collect::<Vec<_>>());
    println!("top: {:?}", top.iter().map(|f| &f.points).collect::<Vec<_>>());

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
    essential
        .on(&right, Ebc::T(|_| 273.15))
        .on(&bot, Ebc::T(|_| 273.15))
        .on(&top, Ebc::T(|_| 273.15));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&left_flux, Nbc::Qt(|_| 5e5));

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
