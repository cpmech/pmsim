#![allow(unused)]

use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::analytical::ElastPlaneStrainPresCylin;
use pmsim::prelude::*;
use pmsim::util::compare_results;
use russell_lab::*;

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

const NAME: &str = "test_spo_751_pres_cylin";
const GENERATE_MESH: bool = false;
const SAVE_FIGURE: bool = true;

const A: f64 = 100.0; // inner radius
const B: f64 = 200.0; // outer radius
const P_ARRAY: [f64; 6] = [0.0, 0.1, 0.14, 0.18, 0.19, 0.192]; // inner pressure
const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_751_press_cylin() -> Result<(), StrError> {
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

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P_ARRAY[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_incremental(P_ARRAY.len());

    // FEM state
    let mut state = FemState::new(&mesh, &base, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, NAME, None)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // compare the results with Ref #1
    let tol_displacement = 1e-12;
    let tol_stress = 1e-14;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        "spo_751_pres_cylin.json",
        tol_displacement,
        tol_stress,
        1,
    )?;
    assert!(all_good);

    // compare with analytical solution
    // let ana = ElastPlaneStrainPresCylin::new(A, B, P_ARRAY[P_ARRAY.len() - 1], 0.0, YOUNG, POISSON)?;
    // let r = mesh.points[ref_point_id].coords[0];
    // assert_eq!(mesh.points[ref_point_id].coords[1], 0.0);
    // let eq = base.equations.eq(ref_point_id, Dof::Ux).unwrap();
    // let numerical_ur = state.uu[eq];
    // let error = f64::abs(numerical_ur - ana.ur(r));
    // println!("\nnumerical_ur = {:?}, error = {:?}", numerical_ur, error);
    // approx_eq(numerical_ur, ana.ur(r), 1.29e-4);

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // load summary and associated files
    let (file_io, mesh, base) = PostProc::read_summary(DEFAULT_OUT_DIR, NAME)?;

    // boundaries
    let features = Features::new(&mesh, false);
    let outer_point = features.search_point_ids(At::XY(B, 0.0), any_x)?[0];

    // loop over time stations
    let mut outer_ur = vec![0.0; file_io.indices.len()];
    let inner_pp: Vec<_> = P_ARRAY.iter().map(|p| *p).collect();
    for index in &file_io.indices {
        // load state
        let state = PostProc::read_state(&file_io, *index)?;
        assert_eq!(file_io.times[*index], state.t);

        // radial displacement
        let outer_eq = base.equations.eq(outer_point, Dof::Ux)?;
        outer_ur[*index] = state.uu[outer_eq];
    }

    // plot
    if SAVE_FIGURE {
        let mut curve1 = Curve::new();
        curve1.draw(&outer_ur, &inner_pp);
        let mut plot = Plot::new();
        plot.add(&curve1)
            .grid_labels_legend("outer $u_r$", "inner $P$")
            .save(&format!("{}/{}.svg", DEFAULT_TEST_DIR, NAME))?;
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
        if SAVE_FIGURE {
            let mut opt = Figure::new();
            opt.with_point_marker = false;
            // opt.point_dots = true;
            opt.point_ids = true;
            opt.cell_ids = true;
            opt.with_cell_att = false;
            opt.figure_size = Some((1000.0, 1000.0));
            mesh.draw(
                Some(opt),
                &format!("{}/mesh_{}_{}.svg", DEFAULT_TEST_DIR, NAME, k_str),
                |_, _| {},
            )
            .unwrap();
        }

        // write mesh
        mesh.write(&format!("{}/{}_{}.msh", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // write VTU
        mesh.write_vtu(&format!("{}/{}_{}.vtu", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read(&format!("data/meshes/{}_{}.msh", NAME, k_str)).unwrap()
    }
}
