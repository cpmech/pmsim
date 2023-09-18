use gemlab::prelude::*;
use plotpy::Curve;
use plotpy::Plot;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{format_nanoseconds, Stopwatch};
use std::env;

const NAME: &str = "pressurized_cylinder2d_elast_";
const SAVE_FIGURE_MESH: bool = false;

// constants
const R1: f64 = 3.0; // inner radius
const R2: f64 = 6.0; // outer radius
const P1: f64 = 200.0; // inner pressure (magnitude)
const P2: f64 = 100.0; // outer pressure (magnitude)
const YOUNG: f64 = 1000.0; // Young's modulus
const POISSON: f64 = 0.25; // Poisson's coefficient

/// Calculates the analytical solution (elastic pressurized cylinder)
/// Reference (page 160)
/// Sadd MH (2005) Elasticity: Theory, Applications and Numerics, Elsevier, 474p
struct AnalyticalSolution {
    aa: f64,
    bb: f64,
    c1: f64,
    c2: f64,
}

impl AnalyticalSolution {
    pub fn new() -> Self {
        let rr1 = R1 * R1;
        let rr2 = R2 * R2;
        let drr = rr2 - rr1;
        let dp = P2 - P1;
        let aa = rr1 * rr2 * dp / drr;
        let bb = (rr1 * P1 - rr2 * P2) / drr;
        let c1 = (1.0 + POISSON) / YOUNG;
        let c2 = 1.0 - 2.0 * POISSON;
        AnalyticalSolution { aa, bb, c1, c2 }
    }
    pub fn radial_displacement(&self, r: f64) -> f64 {
        self.c1 * (r * self.c2 * self.bb - self.aa / r)
    }
}

fn main() -> Result<(), StrError> {
    // arguments
    let args: Vec<String> = env::args().collect();
    let kind = if args.len() >= 2 {
        GeoKind::from(&args[1])?
    } else {
        GeoKind::Tri3
    };
    let str_kind = kind.to_string();

    // sizes
    let sizes = if kind.class() == GeoClass::Tri {
        vec![(2, 4), (5, 10), (20, 40), (50, 100), (120, 220)]
    } else {
        vec![(1, 2), (2, 4), (4, 8), (8, 16), (10, 20), (16, 32), (32, 64), (50, 100)]
    };

    // analytical solution
    let ana = AnalyticalSolution::new();

    // numerical solution arrays
    let n = sizes.len();
    let mut results = ConvergenceResults::new(n);

    // print header
    println!("running with {}", str_kind);
    println!("{:>15} {:>6} {:>8}", "TIME", "NDOF", "ERROR");

    // loop over mesh sizes
    let mut idx = 0;
    for (nr, na) in &sizes {
        // mesh
        let mesh = if kind.class() == GeoClass::Tri {
            let delta_x = (R2 - R1) / (*nr as f64);
            let global_max_area = Some(delta_x * delta_x / 2.0);
            Unstructured::quarter_ring_2d(R1, R2, *nr, *na, kind, global_max_area).unwrap()
        } else {
            Structured::quarter_ring_2d(R1, R2, *nr, *na, kind).unwrap()
        };

        // features
        let feat = Features::new(&mesh, false);
        let bottom = feat.search_edges(At::Y(0.0), any_x)?;
        let left = feat.search_edges(At::X(0.0), any_x)?;
        let inner_circle = feat.search_edges(At::Circle(0.0, 0.0, 3.0), any_x)?;
        let outer_circle = feat.search_edges(At::Circle(0.0, 0.0, 6.0), any_x)?;

        // reference point to compare analytical vs numerical result
        let ref_point_id = feat.search_point_ids(At::XY(R1, 0.0), any_x)?[0];

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
        config.control.verbose_timesteps = false;

        // total number of DOF
        let ndof = data.equations.n_equation;
        let str_ndof = ndof.to_string();

        // save mesh figure
        if SAVE_FIGURE_MESH && idx < 4 {
            let suffix = [str_kind.as_str(), "_", str_ndof.as_str()].concat();
            let fn_svg_mesh = ["/tmp/pmsim/", NAME, suffix.as_str(), "dof.svg"].concat();
            let mut fig = Figure::new();
            let mut curve = Curve::new();
            let mut plot = Plot::new();
            let x = mesh.points[ref_point_id].coords[0];
            let y = mesh.points[ref_point_id].coords[1];
            mesh.draw_cells(&mut fig, true).unwrap();
            let mut with_points = true;
            if kind == GeoKind::Tri3 && idx > 2 {
                with_points = false;
            }
            if kind == GeoKind::Tri6 && idx > 1 {
                with_points = false;
            }
            if with_points {
                mesh.draw_point_dots(&mut fig);
            }
            curve
                .set_line_color("red")
                .set_marker_color("red")
                .set_marker_size(10.0)
                .set_marker_style("o");
            curve.draw(&[x], &[y]);
            plot.add(&curve);
            plot.set_equal_axes(true)
                .set_range(-0.1, R2 + 0.1, -0.1, R2 + 0.1)
                .set_figure_size_points(800.0, 800.0)
                .set_labels("x", "y")
                .save(&fn_svg_mesh)
                .unwrap();
        }

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
        let mut stopwatch = Stopwatch::new("");
        sim.run(&mut state)?;
        results.time[idx] = stopwatch.stop();

        // compute error
        let r = mesh.points[ref_point_id].coords[0];
        let eq = data.equations.eq(ref_point_id, Dof::Ux).unwrap();
        let numerical_ur = state.uu[eq];
        let error = f64::abs(numerical_ur - ana.radial_displacement(r));

        // results
        results.name = str_kind.clone();
        results.ndof[idx] = ndof;
        results.error[idx] = error;
        let ns = format_nanoseconds(results.time[idx]);
        println!("{:>15} {:>6} {:>8.2e}", ns, ndof, error);

        // next mesh
        idx += 1;
    }

    // save results
    let fn_results = ["/tmp/pmsim/", NAME, str_kind.as_str(), "_results.json"].concat();
    results.write(&fn_results)?;
    Ok(())
}
