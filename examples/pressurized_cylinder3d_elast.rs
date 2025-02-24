use gemlab::prelude::*;
use plotpy::Curve;
use plotpy::Surface;
use pmsim::analytical::ElastPlaneStrainPresCylin;
use pmsim::prelude::*;
use pmsim::util::ConvergenceResults;
use russell_lab::*;
use russell_sparse::Genie;
use std::env;

const NAME: &str = "pressurized_cylinder3d_elast";
const SAVE_FIGURE_MESH: bool = false;
const WRITE_K: bool = false;

// constants
const NPOINT_MAX: usize = 1_200_000; // max npoint
const THICKNESS: f64 = 1.0; // out-of-plane thickness
const NZ: usize = 1; // out-of-plane divisions
const R1: f64 = 3.0; // inner radius
const R2: f64 = 6.0; // outer radius
const P1: f64 = 200.0; // inner pressure (magnitude)
const P2: f64 = 100.0; // outer pressure (magnitude)
const YOUNG: f64 = 1000.0; // Young's modulus
const POISSON: f64 = 0.25; // Poisson's coefficient

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
        GeoKind::Tet4
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
            vec![(5, 10), (20, 40), (30, 70)]
        } else {
            // vec![(50, 85)] // very fine mesh for paper
            vec![(2, 4), (10, 20), (20, 40)]
        }
    } else {
        if kind == GeoKind::Hex8 {
            // vec![(88, 100)] // good
            // vec![(89, 100)] // bad
            vec![(4, 8), (12, 16), (40, 50), (120, 180)]
        } else {
            vec![(1, 2), (2, 4), (4, 8), (8, 16), (10, 20), (16, 32), (32, 64), (50, 100)]
        }
    };

    // analytical solution
    let ana = ElastPlaneStrainPresCylin::new(R1, R2, P1, P2, YOUNG, POISSON)?;

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
            let den = if kind == GeoKind::Tet4 { 6.0 } else { 2.0 };
            let global_max_volume = Some(delta_x * delta_x * delta_x / den);
            // println!("0. max vol = {:?}", global_max_volume);
            Unstructured::quarter_ring_3d(R1, R2, THICKNESS, *nr, *na, kind, global_max_volume, true).unwrap()
        } else {
            let wr = vec![1.0; *nr];
            Structured::quarter_ring_3d(R1, R2, THICKNESS, &wr, *na, NZ, kind, true).unwrap()
        };

        // println!("1. npoint = {}, ncell = {}", mesh.points.len(), mesh.cells.len());

        // check mesh
        if mesh.points.len() > NPOINT_MAX {
            return Err("too many points");
        }
        mesh.check_all().unwrap();
        mesh.check_overlapping_points(0.001).unwrap();

        // println!("2. mesh verified");

        // features
        let features = Features::new(&mesh, false);
        let faces_y_min = features.search_faces(At::Y(0.0), any_x)?;
        let faces_x_min = features.search_faces(At::X(0.0), any_x)?;
        let faces_inner = features.search_faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, R1), any_x)?;
        let faces_outer = features.search_faces(At::Cylinder(0.0, 0.0, 0.0, 0.0, 0.0, 1.0, R2), any_x)?;
        let faces_z_min = features.search_faces(At::Z(0.0), any_x)?;
        let faces_z_max = features.search_faces(At::Z(THICKNESS), any_x)?;

        // check boundaries
        if kind == GeoKind::Hex8 {
            assert_eq!(faces_inner.all.len(), *na * NZ);
            assert_eq!(faces_outer.all.len(), *na * NZ);
        }

        // println!("3. found boundaries");

        // reference point to compare analytical vs numerical result
        let ref_point_id = features.search_point_ids(At::XYZ(R1, 0.0, 0.0), any_x)?[0];
        array_approx_eq(&mesh.points[ref_point_id].coords, &[R1, 0.0, 0.0], 1e-15);

        // study point (for debugging)
        let study_point = features.search_point_ids(At::XYZ(0.0, R2, 0.0), any_x)?[0];
        array_approx_eq(&mesh.points[study_point].coords, &[0.0, R2, 0.0], 1e-13); // << some error

        // parameters
        let param1 = ParamSolid {
            density: 1.0,
            stress_strain: StressStrain::LinearElastic {
                young: YOUNG,
                poisson: POISSON,
            },
            ngauss: None,
        };
        let base = FemBase::new(&mesh, [(1, Elem::Solid(param1))])?;

        // total number of DOF
        let ndof = base.dofs.size();
        let n_str = format!("{:0>5}", ndof);

        // println!("4. NDOF = {}", ndof);

        // filepaths
        let ext = if ndof < 20000 { "svg" } else { "png" };
        let path_mesh = format!("/tmp/pmsim/{}_{}_{}_{}.{}", NAME, g_str, k_str, n_str, ext);

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
            cylin_in.set_surf_color("#ff000020");
            cylin_out.set_surf_color("#ff000020");
            cylin_in.draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], R1, 5, 81)?;
            cylin_out.draw_cylinder(&[0.0, 0.0, 0.0], &[0.0, 0.0, 1.0], R2, 5, 81)?;

            // figure settings
            let mut fig = Figure::new();
            fig.size(800.0, 800.0)
                .canvas_points()
                .set_marker_size(2.5)
                .set_marker_line_color("black");
            fig.show_point_dots(if ndof < 3100 { true } else { false })
                .extra(|plot, before| {
                    if !before {
                        plot.add(&cylin_in);
                        plot.add(&cylin_out);
                        plot.add(&curve);
                    }
                })
                .draw(&mesh, &path_mesh)?;
        }

        // essential boundary conditions
        let mut essential = Essential::new();
        essential
            .faces(&faces_x_min, Dof::Ux, 0.0)
            .faces(&faces_y_min, Dof::Uy, 0.0)
            .faces(&faces_z_min, Dof::Uz, 0.0)
            .faces(&faces_z_max, Dof::Uz, 0.0);

        // natural boundary conditions
        let mut natural = Natural::new();
        natural
            .faces(&faces_inner, Nbc::Qn, -P1)
            .faces(&faces_outer, Nbc::Qn, -P2);

        // configuration
        let mut config = Config::new(&mesh);
        config
            .set_linear_problem(true)
            .set_verbose_timesteps(false)
            .set_save_vismatrix_file(false)
            .set_save_matrix_market_file(WRITE_K)
            .set_lin_sol_genie(genie)
            .access_lin_sol_params()
            .umfpack_enforce_unsymmetric_strategy = enforce_unsym_strategy;

        // FEM state
        let mut state = FemState::new(&mesh, &base, &essential, &config)?;

        // File IO
        let mut file_io = FileIo::new();

        // println!("5. running simulation");

        // solution
        let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
        let mut stopwatch = Stopwatch::new();
        match solver.solve(&mut state, &mut file_io) {
            Err(e) => {
                println!("{:?} failed with: {}", genie, e);
                continue;
            }
            Ok(..) => (),
        }
        results.time[idx] = stopwatch.stop();

        // println!("5. computing error");

        // compute error
        let r = mesh.points[ref_point_id].coords[0];
        assert_eq!(mesh.points[ref_point_id].coords[1], 0.0);
        let eq = base.dofs.eq(ref_point_id, Dof::Ux).unwrap();
        let numerical_ur = state.u[eq];
        let error = f64::abs(numerical_ur - ana.ur(r));

        // study point error
        let eq = base.dofs.eq(study_point, Dof::Uy).unwrap();
        let numerical_ur = state.u[eq];
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

        // next mesh
        idx += 1;
    }

    // save results
    results.write_json(&path_json)?;
    Ok(())
}
