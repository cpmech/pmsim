#![allow(unused)]

use gemlab::prelude::*;
use pmsim::analytical::ElastPlaneStrainPresCylin;
use pmsim::prelude::*;
use pmsim::util::verify_results;
use russell_lab::*;

const NAME: &str = "test_spo_751_pres_cylin";
const GENERATE_MESH: bool = true;
const SAVE_FIGURE: bool = true;

const R1: f64 = 100.0; // inner radius
const R2: f64 = 200.0; // outer radius
const P1: [f64; 6] = [0.0, 0.1, 0.14, 0.18, 0.19, 0.192]; // inner pressure
const P2: f64 = 0.0; // outer pressure (magnitude)
const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_751_press_cylin() -> Result<(), StrError> {
    // mesh
    let kind = GeoKind::Qua4;
    let mesh = generate_or_read_mesh(1, kind, GENERATE_MESH);

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, R1), any_x)?;

    // reference point to compare analytical vs numerical result
    let ref_point_id = features.search_point_ids(At::XY(R1, 0.0), any_x)?[0];
    array_approx_eq(&mesh.points[ref_point_id].coords, &[R1, 0.0], 1e-15);

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
    natural.edges_fn(&inner_circle, Nbc::Qn, |t| -P1[t as usize]);

    // configuration
    let mut config = Config::new(&mesh);
    config.set_incremental(P1.len());

    // FEM state
    let mut state = FemState::new(&mesh, &base, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, NAME, None)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // compute error
    let ana = ElastPlaneStrainPresCylin::new(R1, R2, P1[P1.len() - 1], P2, YOUNG, POISSON)?;
    let r = mesh.points[ref_point_id].coords[0];
    assert_eq!(mesh.points[ref_point_id].coords[1], 0.0);
    let eq = base.equations.eq(ref_point_id, Dof::Ux).unwrap();
    let numerical_ur = state.uu[eq];
    let error = f64::abs(numerical_ur - ana.ur(r));
    println!("\nnumerical_ur = {:?}, error = {:?}", numerical_ur, error);
    // approx_eq(numerical_ur, ana.ur(r), 1.29e-4);

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // load essential files
    let (file_io, mesh, base) = PostProc::read_essential(DEFAULT_OUT_DIR, NAME)?;

    // boundaries
    let features = Features::new(&mesh, false);
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, R1), any_x)?;
    // let inner_cell = inner_circle.all[0].

    // loop over time stations
    for index in &file_io.indices {
        // load state
        let state = PostProc::read_state(&file_io, *index)?;
        assert_eq!(file_io.times[*index], state.t);

        // get stress
    }

    Ok(())
}

/// Generate or read mesh
fn generate_or_read_mesh(att: usize, kind: GeoKind, generate: bool) -> Mesh {
    let k_str = kind.to_string();

    if generate {
        // generate mesh
        let wr = &[16.0, 20.0, 28.0, 36.0];
        let na = 3;
        let mesh = Structured::quarter_ring_2d(R1, R2, wr, na, GeoKind::Qua8, true).unwrap();
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
        mesh.write_json(&format!("{}/{}_{}.json", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // write VTU
        mesh.write_vtu(&format!("{}/{}_{}.vtu", DEFAULT_TEST_DIR, NAME, k_str))
            .unwrap();

        // return mesh
        mesh
    } else {
        // read mesh
        Mesh::read_json(&format!("data/meshes/{}_{}.json", NAME, k_str)).unwrap()
    }
}
