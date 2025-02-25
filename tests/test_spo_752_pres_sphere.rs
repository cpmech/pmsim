use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresSphere};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::PI;
use russell_lab::{approx_eq, array_approx_eq};

// This test runs the Example 7.5.2 (aka 752) on page 247 of Ref #1 (aka SPO's book)
//
// Nonetheless, here we use a quarter-ring geometry instead of a 1/12 slice.
//
// This is an axisymmetric problem; hence the ring becomes an octant of a spherical shell.
//
// Axisymmetric
// y ^
//   |
//   ***=---__
//   |        '*._        A slice of an octant of a spherical shell
//   |            *._
//   |               *.
//   ***=-__           *.
//   .      '-.          *
//             *.         *
//   .        P  *         *
//                *         *
//   .             *         *
//                 #         #
//   o -   -   -   # ------- # --> x
//                 a         b
//
// # Reference
//
// 1. de Souza Neto EA, Peric D, Owen DRJ (2008) Computational methods for plasticity,
//    Theory and applications, Wiley, 791p

const NAME_MESH: &str = "spo_751_pres_cylin"; // same as 751
const NAME_COLLAPSE: &str = "spo_752_pres_sphere_collapse";
const NAME_RESIDUAL: &str = "spo_752_pres_sphere_residual";
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius

const P_ARRAY_COLLAPSE: [f64; 5] = [0.0, 0.15, 0.3, 0.33, 0.33269]; // inner pressure
const P_MAX_RES: f64 = 0.28; // maximum pressure achieved by the residual simulation before unloading completely to zero
const P_ARRAY_RESIDUAL: [f64; 4] = [0.0, 0.15, P_MAX_RES, 0.0];
const P_SHOW_COLLAPSE: [f64; 3] = [0.15, 0.30, 0.33];

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 0.24; // uniaxial yield strength (= σy_spo due to axisymmetry)
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_752_pres_sphere() -> Result<(), StrError> {
    // mesh
    let kind = GeoKind::Qua4;
    let mesh = Mesh::read(&format!("data/spo/{}_{}.msh", NAME_MESH, kind.to_string())).unwrap();

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, A), any_x)?;

    // reference point to compare analytical vs numerical result
    let ref_point_id = features.search_point_ids(At::XY(A, 0.0), any_x)?[0];
    array_approx_eq(&mesh.points[ref_point_id].coords, &[A, 0.0], 1e-15);

    // parameters
    let param1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: 0.0,
            z_ini: 0.24,
        },
        ngauss: Some(NGAUSS),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(param1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&left, Dof::Ux, 0.0).edges(&bottom, Dof::Uy, 0.0);

    // run the collapse test
    run_test(false, &mesh, &base, &essential, &inner_circle)?;

    // run the residual stress test
    run_test(true, &mesh, &base, &essential, &inner_circle)?;
    Ok(())
}

fn run_test(
    residual: bool,
    mesh: &Mesh,
    base: &FemBase,
    essential: &Essential,
    inner_circle: &Edges,
) -> Result<(), StrError> {
    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        .set_tol_mdu_rel(1e-7)
        .set_tol_rr_abs(1e-7)
        .set_lagrange_mult_method(false)
        .set_arc_length_method(false)
        .set_ini_trial_load_factor(0.05)
        .update_model_settings(1)
        .set_save_strain(true);

    // natural boundary conditions and configuration
    let mut natural = Natural::new();
    let name = if residual {
        natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY_RESIDUAL[t as usize]);
        config.set_incremental(P_ARRAY_RESIDUAL.len());
        NAME_RESIDUAL
    } else {
        natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY_COLLAPSE[t as usize]);
        config.set_incremental(P_ARRAY_COLLAPSE.len());
        NAME_COLLAPSE
    };

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", name)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // compare the results with Ref #1
    let tol_displacement = if residual { 1.78e-2 } else { 4.33e-2 };
    let tol_stress = if residual { 2.41e-2 } else { 2.55e-2 };
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", name),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    // analyze results
    analyze_results(residual)?;
    Ok(())
}

fn analyze_results(residual: bool) -> Result<(), StrError> {
    // select name and loading array
    let (name, p_array) = if residual {
        (NAME_RESIDUAL, Vec::from(&P_ARRAY_RESIDUAL))
    } else {
        (NAME_COLLAPSE, Vec::from(&P_ARRAY_COLLAPSE))
    };

    // load summary and associated files
    let (post, mut memo) = PostProc::new("/tmp/pmsim", name)?;
    let mesh = post.mesh();
    let base = post.base();

    // boundaries
    let features = Features::new(mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);
    let eq_ux = base.dofs.eq(outer_point, Dof::Ux)?;

    // analytical solution
    let ana = PlastPlaneStrainPresSphere::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut outer_ur = vec![0.0; post.n_state()];
    let inner_pp: Vec<_> = p_array.iter().map(|p| *p).collect();
    let mut first_rr = true;
    let mut rr = Vec::new();
    let mut pp_arr = Vec::new();
    let mut sh_arr = Vec::new();
    let mut sr_arr = Vec::new();
    for index in 0..post.n_state() {
        // load state
        let state = post.read_state(index)?;

        // radial displacement
        let ub_num = state.u[eq_ux];
        outer_ur[index] = ub_num;

        // get stresses
        let res = post.gauss_stresses(&mut memo, &state, &lower_cells, |x, y, _| {
            let alpha = f64::atan2(y, x) * 180.0 / PI;
            alpha < 15.0
        })?;

        // convert to polar coordinates and compare with analytical solution
        let pp = p_array[index];
        if !residual && P_SHOW_COLLAPSE.contains(&pp) || residual && index == post.n_state() - 1 {
            pp_arr.push(pp);
            sh_arr.push(Vec::new());
            sr_arr.push(Vec::new());
            for i in 0..res.xx.len() {
                let (r, sr, sh, _) = cartesian_to_polar(res.xx[i], res.yy[i], res.txx[i], res.tyy[i], res.txy[i]);
                if first_rr {
                    rr.push(r);
                }
                sh_arr.last_mut().unwrap().push(sh);
                sr_arr.last_mut().unwrap().push(sr);
                if residual {
                    let (sr_ana, sh_ana) = ana.calc_sr_sh_residual(r, P_MAX_RES)?;
                    approx_eq(sr, sr_ana, 0.003);
                    approx_eq(sh, sh_ana, 0.003);
                } else {
                    let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                    approx_eq(sr, sr_ana, 0.0036);
                    approx_eq(sh, sh_ana, 0.0077);
                }
            }
            first_rr = false;
        }
    }

    // plot
    if SAVE_FIGURE {
        let mut plot = ana.plot_results(&pp_arr, residual, P_MAX_RES, |plot, index| {
            if index == 0 {
                // load-displacement curve
                let mut curve = Curve::new();
                curve
                    .set_line_style("--")
                    .set_marker_style("o")
                    .set_marker_void(true)
                    .set_label("numerical")
                    .draw(&outer_ur, &inner_pp);
                plot.add(&curve);
            } else if index == 1 {
                // hoop stress-strain curve
                let mut curve = Curve::new();
                curve
                    .set_label("Gauss points")
                    .set_line_style("None")
                    .set_marker_style(".")
                    .set_marker_color("black")
                    .set_marker_line_color("black");
                for i in 0..sh_arr.len() {
                    curve.draw(&rr, &sh_arr[i]);
                }
                plot.add(&curve);
            } else if index == 2 {
                // radial stress-strain curve
                let mut curve = Curve::new();
                curve
                    .set_label("numerical")
                    .set_line_style("None")
                    .set_marker_style(".")
                    .set_marker_color("black")
                    .set_marker_line_color("black");
                for i in 0..sr_arr.len() {
                    curve.draw(&rr, &sr_arr[i]);
                }
                plot.add(&curve);
            } else if index == 3 {
                // legend
                let mut curve = Curve::new();
                curve
                    .set_label("numerical")
                    .set_line_style("None")
                    .set_marker_style(".")
                    .set_marker_color("black")
                    .set_marker_line_color("black")
                    .draw(&[0], &[0]);
                plot.add(&curve);
            }
        });
        plot.set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", name))?;
    }

    Ok(())
}
