use gemlab::prelude::*;
use plotpy::Curve;
use pmsim::analytical::{cartesian_to_polar, PlastPlaneStrainPresSphere};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::array_approx_eq;
use russell_lab::base::read_data;
use russell_lab::math::PI;

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

const NAME: &str = "spo_752_pres_sphere";
const NAME_MESH: &str = "spo_751_pres_cylin"; // same as 751
const SAVE_FIGURE: bool = false;
const VERBOSE_LEVEL: usize = 0;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius
const P_ARRAY: [f64; 5] = [0.0, 0.15, 0.3, 0.33, 0.33269]; // inner pressure
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

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        .set_lagrange_mult_method(false)
        .set_incremental(P_ARRAY.len());

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // compare the results with Ref #1
    let tol_displacement = 2e-1; // TODO: check why the difference is too high
    let tol_stress = 1e-1; // TODO: check why the difference is too high
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // load summary and associated files
    let (file_io, mesh, base) = PostProc::read_summary("/tmp/pmsim", NAME)?;
    let mut post = PostProc::new(&mesh, &base);

    // boundaries
    let features = Features::new(&mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let lower_cells = features.get_cells_via_2d_edges(&bottom);

    // analytical solution
    let ana = PlastPlaneStrainPresSphere::new(A, B, YOUNG, POISSON, Y).unwrap();

    // loop over time stations
    let mut outer_ur = vec![0.0; file_io.indices.len()];
    let inner_pp: Vec<_> = P_ARRAY.iter().map(|p| *p).collect();
    let mut rr = Vec::new();
    let mut sh_p15 = Vec::new();
    let mut sr_p15 = Vec::new();
    let mut sh_p30 = Vec::new();
    let mut sr_p30 = Vec::new();
    for index in &file_io.indices {
        // load state
        let state = PostProc::read_state(&file_io, *index)?;
        assert_eq!(file_io.times[*index], state.t);

        // radial displacement
        let outer_eq = base.equations.eq(outer_point, Dof::Ux)?;
        let ub_num = state.uu[outer_eq];
        outer_ur[*index] = ub_num;

        // get stresses
        let pp = P_ARRAY[*index];
        if pp == 0.15 || pp == 0.3 {
            let res = post.gauss_stresses(&lower_cells, &state, |x, y, _| {
                let alpha = f64::atan2(y, x) * 180.0 / PI;
                alpha < 15.0
            })?;
            for i in 0..res.xx.len() {
                let (r, sr, sh, _) = cartesian_to_polar(res.xx[i], res.yy[i], res.txx[i], res.tyy[i], res.txy[i]);
                if pp == 0.15 {
                    rr.push(r);
                    sh_p15.push(sh);
                    sr_p15.push(sr);
                } else if pp == 0.30 {
                    sh_p30.push(sh);
                    sr_p30.push(sr);
                }
                // let (sr_ana, sh_ana) = ana.calc_sr_sh(r, pp)?;
                // approx_eq(sr, sr_ana, 0.0036);
                // approx_eq(sh, sh_ana, 0.0063);
            }
        }
    }

    // plot
    if SAVE_FIGURE {
        let spo = read_data("data/spo/spo_752_pres_sphere_pp_vs_ub.tsv", &["x", "y"])?;
        let mut plot = ana.plot_results(&[0.15, 0.30], |plot, index| {
            if index == 0 {
                let mut curve_spo = Curve::new();
                curve_spo
                    .set_label("de Souza Neto et al. (2008)")
                    .set_line_style("None")
                    .set_marker_style("D")
                    .set_marker_void(true)
                    .draw(&spo["x"], &spo["y"]);
                plot.add(&curve_spo);
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
                    .draw(&rr, &sh_p30);
                plot.add(&curve);
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("s")
                    .set_marker_void(true)
                    .draw(&rr, &sh_p15);
                plot.add(&curve);
            } else if index == 2 {
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("o")
                    .set_marker_void(true)
                    .draw(&rr, &sr_p30);
                plot.add(&curve);
                let mut curve = Curve::new();
                curve
                    .set_line_style("None")
                    .set_marker_style("s")
                    .set_marker_void(true)
                    .draw(&rr, &sr_p15);
                plot.add(&curve);
            }
        });
        plot.set_figure_size_points(600.0, 450.0)
            .save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}
