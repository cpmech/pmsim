use gemlab::prelude::*;
use plotpy::{Curve, Plot, SlopeIcon};
use pmsim::{prelude::*, StrError};

fn main() -> Result<(), StrError> {
    // constants
    const R1: f64 = 3.0; // inner radius
    const R2: f64 = 6.0; // outer radius
    const P1: f64 = 200.0; // inner pressure (magnitude)
    const P2: f64 = 100.0; // outer pressure (magnitude)
    const YOUNG: f64 = 1000.0;
    const POISSON: f64 = 0.25;

    // let sizes = &[(1, 2), (2, 4), (4, 8), (8, 16), (10, 20), (16, 32), (32, 64), (50, 100)];
    let sizes = &[(2, 4), (5, 10), (20, 40), (50, 100), (120, 220)];
    let n = sizes.len();
    let mut arr_ndof = vec![0.0; n];
    let mut arr_error = vec![0.0; n];

    let mut idx = 0;
    for (nr, na) in sizes {
        // mesh
        // let global_max_area = 2.0 / ((nr * na) as f64);
        let delta_x = (R2 - R1) / (*nr as f64);
        let global_max_area = delta_x * delta_x / 2.0;
        let mesh = Unstructured::quarter_ring_2d(R1, R2, *nr, *na, false, Some(global_max_area), "").unwrap();

        // write msh and svg files
        if false {
            let out_dir = "/tmp/pmsim/";
            let base = "quarter_ring2d_tri3_";
            let pps = mesh.points.len().to_string();
            let ccs = mesh.cells.len().to_string();
            let fn_msh = [out_dir, base, pps.as_str(), "points_", ccs.as_str(), "cells.msh"].concat();
            let fn_svg = [out_dir, base, pps.as_str(), "points_", ccs.as_str(), "cells.svg"].concat();
            let ncell = mesh.cells.len();
            if ncell < 2900 {
                draw_mesh(&mesh, false, ncell < 190, false, &fn_svg)?;
            }
            mesh.write_text_file(&fn_msh)?;
        }

        // features
        let find = Find::new(&mesh, None);
        let bottom = find.edges(At::Y(0.0), any_x)?;
        let left = find.edges(At::X(0.0), any_x)?;
        let inner_circle = find.edges(At::Circle(0.0, 0.0, 3.0), any_x)?;
        let outer_circle = find.edges(At::Circle(0.0, 0.0, 6.0), any_x)?;

        // parameters, DOFs, and configuration
        let param1 = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: YOUNG,
                poisson: POISSON,
            },
        };
        let data = Data::new(&mesh, [(1, Element::Solid(param1))])?;
        let mut config = Config::new();
        config.linear_problem = true;
        config.control.verbose_iterations = false;
        config.control.verbose_timesteps = false;

        // essential boundary conditions
        let mut essential = Essential::new();
        essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0));

        // natural boundary conditions
        let mut natural = Natural::new();
        natural
            .on(&inner_circle, Nbc::Qn(|_| -P1))
            .on(&outer_circle, Nbc::Qn(|_| -P2));

        // simulation state
        let mut state = State::new(&data, &config)?;

        // run simulation
        let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
        sim.run(&mut state)?;

        // analytical solution: radial displacement (ur)
        // Reference (page 160)
        // Sadd MH (2005) Elasticity: Theory, Applications and Numerics, Elsevier, 474p
        let rr1 = R1 * R1;
        let rr2 = R2 * R2;
        let drr = rr2 - rr1;
        let dp = P2 - P1;
        let aa = rr1 * rr2 * dp / drr;
        let bb = (rr1 * P1 - rr2 * P2) / drr;
        let c1 = (1.0 + POISSON) / YOUNG;
        let c2 = 1.0 - 2.0 * POISSON;
        let analytical_ur = |r: f64| c1 * (r * c2 * bb - aa / r);

        let ref_point_id = 0;

        let r = mesh.points[ref_point_id].coords[0];
        let eq = data.equations.eq(ref_point_id, Dof::Ux).unwrap();
        let numerical_ur = state.uu[eq];
        let error = f64::abs(numerical_ur - analytical_ur(r));

        arr_ndof[idx] = data.equations.n_equation as f64;
        arr_error[idx] = error;
        println!("ndof = {:6}, err = {:.2e}", data.equations.n_equation, error);

        idx += 1;
    }

    // let arr_ndof = vec![42.0, 214.0, 2980.0, 18476.0, 105762.0];
    // let arr_error = vec![6.83e-2, 6.91e-3, 5.52e-4, 7.57e-5, 6.17e-6];

    let mut curve = Curve::new();
    let mut icon = SlopeIcon::new();
    icon.set_length(0.265);
    curve.draw(&arr_ndof, &arr_error);
    icon.draw(-1.0, 3.25e2, 0.4e-3);
    let mut plot = Plot::new();
    plot.set_log_x(true)
        .set_log_y(true)
        .add(&curve)
        .add(&icon)
        .grid_and_labels("NDOF", "ERROR");
    plot.save("/tmp/pmsim/pressurized_cylinder2d_elast_tri3_convergence.svg")?;

    Ok(())
}
