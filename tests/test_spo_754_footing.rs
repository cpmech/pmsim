#![allow(unused)]

use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

const SAVE_FIGURE: bool = true;

#[test]
fn test_solid_smith_5d2_tri3_plane_strain() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::from_text_file("data/meshes/test_spo_754_footing.msh")?;
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
        mesh.draw(Some(opt), "/tmp/pmsim/test_spo_754_footing.svg", |_, _| {})?;
    }

    // features
    let feat = Features::new(&mesh, false);
    let left = feat.search_edges(At::X(0.0), any_x)?;
    let right = feat.search_edges(At::X(5.0), any_x)?;
    let bottom = feat.search_edges(At::Y(0.0), any_x)?;
    let top = feat.search_edges(At::Y(5.0), any_x)?;
    let footing = feat.search_edges(At::Y(5.0), |x| x[0] <= 0.5)?;
    // for f in footing {
    //     println!("{:?}", f.points);
    // }

    // E   = 1e+07  // kPa
    // nu  = 0.48   // -
    // qy0 = 848.7  // kPa
    // H   = 0      // kPa
    // rho = 2      // Mg/m3

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: 1e+07,
            poisson: 0.48,
        },
    };
    let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))])?;

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

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .on(&left, Ebc::Ux(|_| 0.0))
        .on(&right, Ebc::Ux(|_| 0.0))
        .on(&bottom, Ebc::Uy(|_| 0.0))
        .on(&footing, Ebc::Uy(|t: f64| UY[t as usize]));

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_tol_rr(1e-6)
        .set_ngauss(att, 4)
        .set_dt(|_| 1.0)
        .set_dt_out(|_| 1.0)
        .set_t_fin((UY.len() - 1) as f64)
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&input, &config)?;
    let mut output = FemOutput::new(&input, None, None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    Ok(())
}
