use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresSphere};
use pmsim::material::{Plotter, PlotterData};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::array_approx_eq;
use russell_lab::base::read_data;
use russell_lab::math::{PI, SQRT_2_BY_3};

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
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 1;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius
const P_ARRAY_COLLAPSE: [f64; 5] = [0.0, 0.15, 0.3, 0.33, 0.33269]; // inner pressure
const P_ARRAY_RESIDUAL: [f64; 4] = [0.0, 0.15, 0.28, 0.0]; // inner pressure
const PA_COLLAPSE: f64 = 0.15;
const PB_COLLAPSE: f64 = 0.30;
const PA_RESIDUAL: f64 = 0.15;
const PB_RESIDUAL: f64 = 0.28;
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
        // stress_strain: StressStrain::LinearElastic {
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

    // run the collapse case
    run_spo_752(&mesh, &inner_circle, &base, &essential, true)?;

    // run the residual case
    // run_spo_752(&mesh, &inner_circle, &base, &essential, false)?;

    Ok(())
}

fn run_spo_752(
    mesh: &Mesh,
    inner_circle: &Edges,
    base: &FemBase,
    essential: &Essential,
    collapse: bool,
) -> Result<(), StrError> {
    // natural boundary conditions
    let mut natural = Natural::new();
    if collapse {
        natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY_COLLAPSE[t as usize]);
    } else {
        natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY_RESIDUAL[t as usize]);
    }

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        // .set_n_max_iterations(10)
        .set_constant_tangent(false)
        .set_lagrange_mult_method(true)
        .update_model_settings(1)
        .set_save_strain(true);
    if collapse {
        config.set_incremental(P_ARRAY_COLLAPSE.len());
    } else {
        config.set_incremental(P_ARRAY_RESIDUAL.len());
    }

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    let name = if collapse {
        "spo_752_pres_sphere_collapse"
    } else {
        "spo_752_pres_sphere_resid_stress"
    };
    file_io.activate(&mesh, &base, "/tmp/pmsim", name)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // compare the results with Ref #1
    let fn_ref = if collapse {
        "data/spo/spo_752_pres_sphere_collapse_ref.json"
    } else {
        "data/spo/spo_752_pres_sphere_resid_stress_ref.json"
    };
    let tol_displacement = 4.33e-2;
    let tol_stress = 2.55e-2;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SPO,
        fn_ref,
        // "data/spo/spo_752_pres_sphere_resid_stress_ref.json",
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    // post-processing
    // post_processing(collapse, name)
    Ok(())
}

fn post_processing(collapse: bool, name: &str) -> Result<(), StrError> {
    // load summary and associated files
    let (file_io, mesh, base) = PostProc::deprecated_read_summary("/tmp/pmsim", name)?;
    let mut post = PostProc::deprecated_new(&mesh, &base);

    // boundaries
    let features = Features::new(&mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);

    // analytical solution
    let ana = PlastPlaneStrainPresSphere::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut outer_ur = vec![0.0; file_io.indices.len()];
    let inner_pp: Vec<_> = if collapse {
        P_ARRAY_COLLAPSE.iter().map(|p| *p).collect()
    } else {
        P_ARRAY_RESIDUAL.iter().map(|p| *p).collect()
    };
    let mut rr = Vec::new();
    let mut sh_ppa = Vec::new();
    let mut sr_ppa = Vec::new();
    let mut sh_ppb = Vec::new();
    let mut sr_ppb = Vec::new();
    for index in &file_io.indices {
        // load state
        let state = PostProc::deprecated_read_state(&file_io, *index)?;
        assert_eq!(file_io.times[*index], state.t);

        // radial displacement
        let outer_eq = base.dofs.eq(outer_point, Dof::Ux)?;
        let ub_num = state.u[outer_eq];
        outer_ur[*index] = ub_num;

        // get stresses
        let pp = if collapse {
            P_ARRAY_COLLAPSE[*index]
        } else {
            P_ARRAY_RESIDUAL[*index]
        };
        if pp == PA_COLLAPSE || pp == PB_COLLAPSE || pp == PA_RESIDUAL || pp == PB_RESIDUAL {
            let res = post.gauss_stresses(&lower_cells, &state, |x, y, _| {
                let alpha = f64::atan2(y, x) * 180.0 / PI;
                alpha < 15.0
            })?;
            for i in 0..res.xx.len() {
                let (r, sr, sh, _) = cartesian_to_polar(res.xx[i], res.yy[i], res.txx[i], res.tyy[i], res.txy[i]);
                if pp == PA_COLLAPSE || pp == PA_RESIDUAL {
                    rr.push(r);
                    sh_ppa.push(sh);
                    sr_ppa.push(sr);
                } else if pp == PB_COLLAPSE || pp == PB_RESIDUAL {
                    sh_ppb.push(sh);
                    sr_ppb.push(sr);
                }
                // let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                // approx_eq(sr, sr_ana, 0.0036);
                // approx_eq(sh, sh_ana, 0.0063);
            }
        }
    }

    // plot
    if SAVE_FIGURE {
        let pa = if collapse { PA_COLLAPSE } else { PA_RESIDUAL };
        let pb = if collapse { PB_COLLAPSE } else { PB_RESIDUAL };
        let spo = read_data("data/spo/spo_752_pres_sphere_pp_vs_ub.tsv", &["ub", "P"])?;
        let mut plot = ana.plot_results(&[pa, pb], |plot, index| {
            if index == 0 {
                if collapse {
                    let mut curve_spo = Curve::new();
                    curve_spo
                        .set_label("de Souza Neto et al. (2008)")
                        .set_line_style("None")
                        .set_marker_style("D")
                        .set_marker_void(true)
                        .draw(&spo["ub"], &spo["P"]);
                    plot.add(&curve_spo);
                }
                let mut curve = Curve::new();
                curve
                    .set_line_style("--")
                    .set_marker_style("o")
                    .set_marker_void(true)
                    .set_label("numerical")
                    .draw(&outer_ur, &inner_pp);
                plot.add(&curve);
            } else if index == 1 {
                let mut curve = Curve::new();
                curve
                    .set_label("Gauss points")
                    .set_line_style("None")
                    .set_marker_style("o")
                    .set_marker_void(true)
                    .draw(&rr, &sh_ppb);
                plot.add(&curve);
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("s")
                    .set_marker_void(true)
                    .draw(&rr, &sh_ppa);
                plot.add(&curve);
            } else if index == 2 {
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("o")
                    .set_marker_void(true)
                    .draw(&rr, &sr_ppb);
                plot.add(&curve);
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("s")
                    .set_marker_void(true)
                    .draw(&rr, &sr_ppa);
                plot.add(&curve);
            }
        });
        plot.set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", name))?;
    }

    Ok(())
}

// #[test]
fn _test_spo_752_pres_sphere_debug() -> Result<(), StrError> {
    // read summary and associated files
    let name = "spo_752_pres_sphere_resid_stress";
    let (file_io, _, _) = PostProc::deprecated_read_summary("/tmp/pmsim", name)?;

    // loop over time stations
    let cell_id = 0;
    let gauss_id = 1;
    let mut local_states = Vec::new();
    for index in &file_io.indices {
        let state = PostProc::deprecated_read_state(&file_io, *index)?;
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
    let p = local_states.len() - 1;
    // let radius_0 = local_states[0].int_vars[0] * SQRT_2_BY_3;
    // let radius_1 = local_states[p].int_vars[0] * SQRT_2_BY_3;
    // plotter.set_oct_circle(radius_0, |_| {});
    // plotter.set_oct_circle(radius_1, |canvas| {
    //     canvas.set_line_style("-");
    // });
    plotter.save("/tmp/pmsim/spo_752_pres_sphere_resid_stress_local_state.svg")?;
    Ok(())
}
