#![allow(unused)]

use gemlab::prelude::*;
use math::{sign, PI};
use plotpy::{linspace, Curve, Plot};
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

const NAME: &str = "test_durand_farias_example4";
const GENERATE_MESH: bool = true;
const SAVE_FIGURE: bool = true;

// constants
const B: f64 = 2.5; // half-width of the flexible footing
const W: f64 = 15.0 * B; // width of the domain
const H: f64 = 15.0 * B; // height of the domain
const QN: f64 = 200.0; // magnitude of the distributed loading
const E: f64 = 1000.0; // Young's modulus
const NU: f64 = 0.25; // Poisson's coefficient

#[test]
fn test_durand_farias_example4() -> Result<(), StrError> {
    // mesh
    let att = 1;
    let kind = GeoKind::Qua4;
    let mesh = generate_or_read_mesh(att, kind, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(W), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let footing = features.search_edges(At::Y(H), |x| x[0] <= B)?;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic { young: E, poisson: NU },
    };
    let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(|_| 0.0))
        .on(&right, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&footing, Nbc::Qn(|_| -QN));

    // configuration
    let mut config = Config::new(&mesh);
    config.set_ngauss(att, 4);

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None, None)?;

    // solution
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // verification
    let f_aux = |theta: f64| QN * (theta + 0.5 * f64::sin(2.0 * theta)) / PI;
    let f_sigma_v = |x: f64, y: f64| {
        assert!(x >= 0.0 && y >= 0.0 && y <= H);
        let z = H - y;
        let d1 = x - B;
        let d2 = B + x;
        let s = sign(d1);
        let (theta_1, theta_2) = if z == 0.0 {
            (s * PI / 2.0, PI / 2.0)
        } else {
            (s * f64::atan(f64::abs(d1) / z), f64::atan(d2 / z))
        };
        let sigma_v = f_aux(theta_1) - f_aux(theta_2); // compressive is negative
        sigma_v
    };

    println!(">>> {}", f_sigma_v(0.0, 9.0 * B));

    let k_str = kind.to_string();
    let mut plot = Plot::new();

    let mut curve1 = Curve::new();
    let ll = linspace(0.0, 15.0, 101);
    let nn: Vec<_> = ll.iter().map(|l| -f_sigma_v(0.0, l * B) / QN).collect();
    curve1.draw(&nn, &ll);

    plot.set_gaps(0.25, 0.0)
        .set_subplot(1, 2, 1)
        .set_title("Durand and Farias Fig 13")
        .add(&curve1)
        .grid_labels_legend("Normalized stress: $-\\sigma_v/q_n$", "Normalized length: $y/B$");

    plot.set_subplot(1, 2, 2)
        .set_title("Stress near the surface with y = H - d");
    let xx = linspace(0.0, 2.0 * B, 101);
    for d in [0.0, 0.1, B / 2.0, B] {
        let mut curve2 = Curve::new();
        curve2.set_label(&format!("d = {}", d));
        let sv: Vec<_> = xx.iter().map(|x| f_sigma_v(*x, H - d)).collect();
        curve2.draw(&xx, &sv);
        plot.add(&curve2);
    }

    plot.grid_labels_legend("$x$", "$\\sigma_v$")
        .set_figure_size_points(800.0, 300.0)
        .save(&format!("{}/{}_{}.svg", DEFAULT_TEST_DIR, NAME, k_str))?;

    // calculate σz along left vertical boundary
    for edge in &footing {
        let (a, b) = (edge.points[0], edge.points[1]);
        let cell = features.all_2d_edges.get(&(a, b)).unwrap();
        let cell_id = cell[0].0; // the first 0 corresponds to the only boundary cell, the second 0 is the first tuple's item
    }

    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(att: usize, kind: GeoKind, generate: bool) -> Mesh {
    let k_str = kind.to_string();

    if generate {
        // generate mesh
        let (xa, xb, xc, na, nb) = (0.0, B, W, 4, 3);
        let (ya, yb, yc, ma, mb) = (0.0, H - B / 2.0, H, 6, 1);
        let mesh = Structured::rectangle(
            xa,
            Some(xb),
            xc,
            na,
            nb,
            &[ya, yb, yc],
            &[ma, mb],
            &[att, att],
            kind,
            true,
        )
        .unwrap();
        mesh.check_all().unwrap();

        // draw figure
        if SAVE_FIGURE {
            let mut opt = Figure::new();
            opt.with_point_marker = false;
            // opt.point_dots = true;
            opt.point_ids = true;
            opt.cell_ids = true;
            opt.figure_size = Some((1000.0, 1000.0));
            mesh.draw(
                Some(opt),
                &format!("{}/mesh_{}_{}.svg", DEFAULT_TEST_DIR, NAME, k_str),
                |_, _| {},
            )
            .unwrap();
        }

        // write mesh
        mesh.write_json(&format!("{}/{}_{}.json", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read_json(&format!("data/meshes/{}_{}.json", NAME, k_str)).unwrap()
    }
}
