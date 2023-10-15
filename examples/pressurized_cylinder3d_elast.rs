use gemlab::prelude::*;
use plotpy::Curve;
use plotpy::Surface;
use pmsim::prelude::*;
use russell_lab::*;
use russell_sparse::Genie;
use std::env;

const NAME: &str = "pressurized_cylinder3d_elast";
const SAVE_FIGURE_MESH: bool = false;
const SAVE_VTU: bool = false;

// constants
const THICKNESS: f64 = 1.0; // out-of-plane thickness
const NZ: usize = 1; // out-of-plane divisions
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
    let genie = if args.len() > 1 {
        Genie::from(&args[1])
    } else {
        Genie::Mumps
    };
    let kind = if args.len() > 2 {
        GeoKind::from(&args[2])?
    } else {
        GeoKind::Hex8
    };
    let enforce_unsym_strategy = if args.len() > 3 {
        args[3].to_lowercase() == "true"
    } else {
        false
    };

    // filepaths
    let g_str = genie.to_string();
    let k_str = kind.to_string();
    let path_json = format!("/tmp/pmsim/{}_{}_{}.json", NAME, g_str, k_str);

    // sizes
    let sizes = if kind.class() == GeoClass::Tet {
        if kind == GeoKind::Tet4 {
            vec![(5, 10), (20, 40), (50, 100)]
        } else {
            vec![(2, 4), (5, 10), (20, 40), (50, 100)]
        }
    } else {
        if kind == GeoKind::Hex8 {
            vec![(4, 8), (12, 16), (40, 50), (120, 180)]
        } else {
            vec![(1, 2), (2, 4), (4, 8), (8, 16), (10, 20), (16, 32), (32, 64), (50, 100)]
        }
    };

    // analytical solution
    let ana = AnalyticalSolution::new();

    // numerical solution arrays
    let n = sizes.len();
    let mut results = ConvergenceResults::new(n);

    // print header
    println!(
        "... running with {:?} ... {:?} ... enforce_unsym_strategy = {} ...",
        kind, genie, enforce_unsym_strategy
    );
    println!(
        "{:>15} {:>6} {:>11} {:>9} {:>10}",
        "TIME", "NDOF", "log10(NDOF)", "ERROR", "UR_STUDY"
    );

    // loop over mesh sizes
    let mut idx = 0;
    for (nr, na) in &sizes {
        let mesh = if kind.class() == GeoClass::Tet {
            let delta_x = (R2 - R1) / (*nr as f64);
            let global_max_volume = Some(delta_x * delta_x * delta_x / 6.0);
            // println!("max vol = {:?}", global_max_volume);
            Unstructured::quarter_ring_3d(R1, R2, THICKNESS, *nr, *na, kind, global_max_volume).unwrap()
        } else {
            Structured::quarter_ring_3d(R1, R2, THICKNESS, *nr, *na, NZ, kind).unwrap()
        };

        // check mesh
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.001).unwrap();

        // features
        let feat = Features::new(&mesh, false);
        let faces_y_min = feat.search_faces(At::Y(0.0), any_x)?;
        let faces_x_min = feat.search_faces(At::X(0.0), any_x)?;
        let faces_inner = feat.search_faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, R1), any_x)?;
        let faces_outer = feat.search_faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, R2), any_x)?;
        let faces_z_min = feat.search_faces(At::Z(0.0), any_x)?;
        let faces_z_max = feat.search_faces(At::Z(THICKNESS), any_x)?;

        // check boundaries
        if kind == GeoKind::Hex8 {
            assert_eq!(faces_inner.len(), *na * NZ);
            assert_eq!(faces_outer.len(), *na * NZ);
        }

        // reference point to compare analytical vs numerical result
        let ref_point_id = feat.search_point_ids(At::XYZ(R1, 0.0, 0.0), any_x)?[0];
        vec_approx_eq(&mesh.points[ref_point_id].coords, &[R1, 0.0, 0.0], 1e-15);

        // study point (for debugging)
        let study_point = feat.search_point_ids(At::XYZ(0.0, R2, 0.0), any_x)?[0];
        vec_approx_eq(&mesh.points[study_point].coords, &[0.0, R2, 0.0], 1e-13); // << some error

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
        config.control.save_vismatrix_file = false;
        config.control.save_matrix_market_file = false;
        config.lin_sol_genie = genie;
        config.lin_sol_params.umfpack_enforce_unsymmetric_strategy = enforce_unsym_strategy;

        // total number of DOF
        let ndof = data.equations.n_equation;
        let n_str = format!("{:0>5}", ndof);

        // filepaths
        let ext = if ndof < 20000 { "svg" } else { "png" };
        let path_mesh = format!("/tmp/pmsim/{}_{}_{}_{}.{}", NAME, g_str, k_str, n_str, ext);
        let path_vtu = format!("/tmp/pmsim/{}_{}_{}_{}.vtu", NAME, g_str, k_str, n_str);

        // save mesh figure
        if SAVE_FIGURE_MESH {
            // reference point
            let mut curve = Curve::new();
            let x = mesh.points[ref_point_id].coords[0];
            let y = mesh.points[ref_point_id].coords[1];
            curve
                .set_line_color("red")
                .set_marker_color("None")
                .set_marker_size(10.0)
                .set_marker_style("o");
            curve.draw(&[x], &[y]);

            // cylinders
            let mut cylin_in = Surface::new();
            let mut cylin_out = Surface::new();
            cylin_in.set_solid_color("#ff000020");
            cylin_out.set_solid_color("#ff000020");
            cylin_in.draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], R1, 5, 81)?;
            cylin_out.draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], R2, 5, 81)?;

            // figure settings
            let mut fig = Figure::new();
            fig.point_ids = false;
            fig.canvas_points.set_marker_size(2.5).set_marker_line_color("black");
            fig.figure_size = Some((800.0, 800.0));
            fig.point_dots = if ndof < 3100 { true } else { false };

            // draw figure
            mesh.draw(Some(fig), &path_mesh, |plot, before| {
                if !before {
                    plot.add(&cylin_in);
                    plot.add(&cylin_out);
                    plot.add(&curve);
                }
            })?;
        }

        // essential boundary conditions
        let mut essential = Essential::new();
        essential
            .on(&faces_x_min, Ebc::Ux(|_| 0.0))
            .on(&faces_y_min, Ebc::Uy(|_| 0.0))
            .on(&faces_z_min, Ebc::Uz(|_| 0.0))
            .on(&faces_z_max, Ebc::Uz(|_| 0.0));

        // natural boundary conditions
        let mut natural = Natural::new();
        natural
            .on(&faces_inner, Nbc::Qn(|_| -P1))
            .on(&faces_outer, Nbc::Qn(|_| -P2));

        // simulation state
        let mut state = State::new(&data, &config)?;

        // run simulation
        let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
        let mut stopwatch = Stopwatch::new("");
        sim.run(&mut state)?;
        results.time[idx] = stopwatch.stop();

        // compute error
        let r = mesh.points[ref_point_id].coords[0];
        assert_eq!(mesh.points[ref_point_id].coords[1], 0.0);
        let eq = data.equations.eq(ref_point_id, Dof::Ux).unwrap();
        let numerical_ur = state.uu[eq];
        let error = f64::abs(numerical_ur - ana.radial_displacement(r));

        // study point error
        let eq = data.equations.eq(study_point, Dof::Uy).unwrap();
        let numerical_ur = state.uu[eq];
        let study_error = numerical_ur; // should be zero with R2 = 2*R1 and P1 = 2*P2

        // results
        results.name = kind.to_string();
        results.ndof[idx] = ndof;
        results.error[idx] = error;
        let ns = format_nanoseconds(results.time[idx]);
        let lx = f64::log10(ndof as f64);
        println!(
            "{:>15} {:>6} {:>11.2} {:>9.2e} {:>10.2e}",
            ns, ndof, lx, error, study_error
        );

        // vtu file
        if SAVE_VTU {
            let post = PostProc::new(&mesh, &feat, &data, &state);
            post.write_vtu(&path_vtu)?;
        }

        // next mesh
        idx += 1;
    }

    // save results
    results.write(&path_json)?;
    Ok(())
}
