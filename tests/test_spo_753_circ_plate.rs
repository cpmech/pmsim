use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::analytical::PlastCircularPlateAxisym;
use pmsim::prelude::*;
use pmsim::util::compare_results;
use russell_lab::*;

const NAME: &str = "spo_753_circ_plate";
const DRAW_MESH_AND_EXIT: bool = false;
const SAVE_FIGURE: bool = true;
const VERBOSE_LEVEL: usize = 0;

const P_ARRAY: [f64; 13] = [
    0.0, 100.0, 200.0, 220.0, 230.0, 240.0, 250.0, 255.0, 257.0, 259.0, 259.5, 259.75, 259.77,
];
const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.24; // Poisson's coefficient
const Z_INI: f64 = 16000.0; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_753_circ_plate() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;
    if DRAW_MESH_AND_EXIT {
        mesh.check_all()?;
        let mut fig = Figure::new();
        return fig
            .show_point_ids(true)
            .size(600.0, 100.0)
            .range_2d(-0.5, 10.5, -0.2, 1.2)
            .draw(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME));
    }

    // features
    let features = Features::new(&mesh, false);
    let top = features.search_edges(At::Y(1.0), any_x)?;
    let right_corner = features.search_point_ids(At::XY(10.0, 0.0), any_x)?;
    assert_eq!(right_corner.len(), 1);

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.point(right_corner[0], Dof::Uy, 0.0);
    // println!("{}", essential);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.edges_fn(&top, Nbc::Qn, |t| -P_ARRAY[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_axisymmetric()
        .set_incremental(P_ARRAY.len())
        .set_lagrange_mult_method(false)
        .set_tol_rr(1e-6)
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // post-processing
    post_processing()?;

    // verify the results
    /*
    let tol_displacement = 1e-1;
    let tol_stress = 1e-1;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);
    */
    Ok(())
}

fn post_processing() -> Result<(), StrError> {
    // load summary and associated files
    let (file_io, mesh, base) = PostProc::read_summary("/tmp/pmsim", NAME)?;
    let post = PostProc::new(&mesh, &base);

    // boundaries
    let features = Features::new(&mesh, false);
    let center = features.search_point_ids(At::XY(0.0, 0.0), any_x)?[0];
    let eq_uy = base.equations.eq(center, Dof::Uy)?;

    // analytical solution
    let ana = PlastCircularPlateAxisym::new(10.0, 1.0, Z_INI);

    // load results
    let mut deflection = vec![0.0; file_io.indices.len()];
    let load: Vec<_> = P_ARRAY.iter().map(|p| *p).collect();
    for index in &file_io.indices {
        let state = PostProc::read_state(&file_io, *index)?;
        deflection[*index] = -state.uu[eq_uy];
    }

    // plot
    if SAVE_FIGURE {
        let ref_xy = Matrix::from_text_file("data/spo/spo_753_plate_deflection_load.tsv")?;
        let ref_x = ref_xy.extract_column(0);
        let ref_y = ref_xy.extract_column(1);
        let mut plot = Plot::new();
        let mut curve_ref = Curve::new();
        curve_ref
            .set_label("de Souza Neto et al.")
            .set_line_style("None")
            .set_marker_style("D")
            .set_marker_void(true)
            .draw(&ref_x, &ref_y);
        let mut curve_num = Curve::new();
        curve_num
            .set_line_style("--")
            .set_marker_style("o")
            .set_marker_void(true)
            .set_label("numerical")
            .draw(&deflection, &load);
        plot.set_horiz_line(ana.get_pp_lim(), "green", ":", 1.0)
            .add(&curve_num)
            .add(&curve_ref)
            .grid_labels_legend("Central deflection", "Distributed load intensity")
            .save(&format!("/tmp/pmsim/{}.svg", NAME))?;
    }

    Ok(())
}
