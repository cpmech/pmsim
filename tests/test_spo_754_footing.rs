#![allow(unused)]

use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, fem::FileIoSummary, prelude::*, util::verify_results};
use russell_lab::*;

const NAME: &str = "test_spo_754_footing";
const SAVE_FIGURE: bool = true;

#[test]
fn test_spo_754_footing() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::from_text_file(&format!("data/meshes/{}.msh", NAME))?;
    let att = mesh.cells[0].attribute;
    if SAVE_FIGURE {
        mesh.check_all()?;
        let mut opt = Figure::new();
        opt.with_point_marker = false;
        // opt.point_dots = true;
        // opt.point_ids = true;
        // opt.cell_ids = true;
        // opt.figure_size = Some((2000.0, 2000.0));
        opt.figure_size = Some((800.0, 800.0));
        mesh.draw(Some(opt), &format!("/tmp/pmsim/{}.svg", NAME), |_, _| {})?;
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(5.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(5.0), any_x)?;
    let footing = features.search_edges(At::Y(5.0), |x| x[0] <= 0.5)?;
    // for f in footing {
    //     println!("{:?}", f.points);
    // }

    // E   = 1e+07  // kPa
    // nu  = 0.48   // -
    // qy0 = 848.7  // kPa
    // H   = 0      // kPa
    // rho = 2      // Mg/m3

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: 1e+07,
            poisson: 0.48,
            z_ini: 848.7,
            hh: 0.0,
        },
    };
    let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))])?;

    /*
    const UY: [f64; 15] = [
        0.0,     //  0
        0.0001,  //  1
        0.00015, //  2
        0.0002,  //  3
        0.00025, //  4
        0.00035, //  5
        0.00045, //  6
        0.00055, //  7
        0.00065, //  8
        0.00075, //  9
        0.00080, // 10
        0.00090, // 11
        0.00110, // 12
        0.00140, // 13
        0.002,   // 14
    ];
    */
    const UY: [f64; 2] = [0.0, 0.0001];

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&footing, Dof::Uy, |t: f64| UY[t as usize]);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_tol_rr(1e-6)
        .set_ngauss(att, 4)
        .set_incremental(UY.len())
        .set_ignore_jacobian_symmetry(true)
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&fem, &config)?;

    // File IO
    let mut file_io = FileIo::new(&fem, Some(NAME.to_string()), None)?;

    // solution
    let mut solver = SolverImplicit::new(&fem, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // analysis
    // let summary = FileIoSummary::read_json(&format!("{}/{}-summary.json", DEFAULT_OUT_DIR, NAME))?;

    // verify the results
    let tol_displacement = 1e-1;
    let tol_stress = 1e+3;
    verify_results(&mesh, NAME, "spo_754_footing.json", tol_displacement, tol_stress, true)?;

    Ok(())
}
