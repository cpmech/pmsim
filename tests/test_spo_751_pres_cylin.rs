use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresCylin};
use pmsim::material::{Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::{PI, SQRT_3};
use russell_lab::{approx_eq, array_approx_eq};

// This test runs the Example 7.5.1 (aka 751) on page 244 of Ref #1 (aka SPO's book)
//
// Nonetheless, here we use a quarter-ring geometry instead of a 1/12 slice.
//
// y ^
//   |
//   ***=---__
//   |        '*._
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

const NAME_MESH: &str = "spo_751_pres_cylin";
const NAME_COLLAPSE: &str = "spo_751_pres_cylin_collapse";
const NAME_RESIDUAL: &str = "spo_751_pres_cylin_residual";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius

const P_ARRAY_COLLAPSE: [f64; 6] = [0.0, 0.1, 0.14, 0.18, 0.19, 0.192]; // inner pressure
const P_MAX_RES: f64 = 0.18; // maximum pressure achieved by the residual simulation before unloading completely to zero
const P_ARRAY_RESIDUAL: [f64; 5] = [0.0, 0.1, 0.14, P_MAX_RES, 0.0];
const P_SHOW_COLLAPSE: [f64; 3] = [0.1, 0.18, 0.19];

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const Y: f64 = 2.0 * 0.24 / SQRT_3; // uniaxial yield strength (2 σy_spo / sq3)
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_751_pres_cylin() -> Result<(), StrError> {
    // mesh
    let kind = GeoKind::Qua4;
    let mesh = generate_or_read_mesh(kind, GENERATE_MESH);

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
        .set_tol_mdu_rel(1e-10)
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
    let tol_displacement = if residual { 1e-13 } else { 1e-11 };
    let tol_stress = 1e-14;
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
    let ana = PlastPlaneStrainPresCylin::new(A, B, YOUNG, POISSON, Y).unwrap();

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

/// Generate or read mesh
fn generate_or_read_mesh(kind: GeoKind, generate: bool) -> Mesh {
    let k_str = kind.to_string();

    if generate {
        // generate mesh
        let wr = &[16.0, 20.0, 28.0, 36.0];
        let na = 3;
        let mesh = Structured::quarter_ring_2d(A, B, wr, na, GeoKind::Qua8, true).unwrap();
        mesh.check_all().unwrap();

        // draw figure
        let mut fig = Figure::new();
        fig.show_point_ids(true)
            .show_cell_ids(true)
            .size(600.0, 600.0)
            .draw(&mesh, &format!("/tmp/pmsim/{}_{}.svg", NAME_MESH, k_str))
            .unwrap();

        // write mesh
        mesh.write(&format!("/tmp/pmsim/{}_{}.msh", NAME_MESH, k_str)).unwrap();

        // write VTU
        mesh.write_vtu(&format!("/tmp/pmsim/{}_{}.vtu", NAME_MESH, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/spo/{}_{}.msh", NAME_MESH, k_str)).unwrap()
    }
}

// #[test]
fn _test_spo_751_pres_cylin_debug() -> Result<(), StrError> {
    // read summary and associated files
    let name = "spo_751_pres_cylin_resid_stress";
    let (post, _) = PostProc::new("/tmp/pmsim", name)?;

    // loop over time stations
    let cell_id = 0;
    let gauss_id = 1;
    let mut local_states = Vec::new();
    for index in 0..post.n_state() {
        let state = post.read_state(index)?;
        println!(
            "t = {:?}",
            state.gauss[cell_id].stress(gauss_id).unwrap().vector().as_data()
        );
        local_states.push(state.gauss[cell_id].get_local_state(gauss_id)?.clone());
    }

    let data = PlotterData::from_states(&local_states);
    let mut plotter = Plotter::new();
    plotter
        .add_2x2(&data, false, |curve, _, _| {
            curve.set_marker_style("o");
        })
        .unwrap();
    let _p = local_states.len() - 1;
    // let radius_0 = local_states[0].int_vars[0] * SQRT_2_BY_3;
    // let radius_1 = local_states[p].int_vars[0] * SQRT_2_BY_3;
    // plotter.set_oct_circle(radius_0, |_| {});
    // plotter.set_oct_circle(radius_1, |canvas| {
    //     canvas.set_line_style("-");
    // });
    plotter.save("/tmp/pmsim/spo_751_pres_cylin_resid_stress_local_state.svg")?;
    Ok(())
}
