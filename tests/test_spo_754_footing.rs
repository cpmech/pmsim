#![allow(unused)]

use gemlab::prelude::*;
use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::util::compare_results;
use russell_lab::*;

const NAME: &str = "test_spo_754_footing";
const SAVE_FIGURE: bool = false;

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_754_footing() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/meshes/{}.msh", NAME))?;
    if SAVE_FIGURE {
        mesh.check_all()?;
        let mut opt = Figure::new();
        opt.with_point_marker = false;
        // opt.point_dots = true;
        opt.point_ids = true;
        opt.cell_ids = true;
        opt.figure_size = Some((1000.0, 1000.0));
        mesh.draw(Some(opt), &format!("/tmp/pmsim/{}.svg", NAME), |plot, _| {
            plot.set_range(0.0, 0.5, 4.5, 5.0);
        })?;
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(500.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let top = features.search_edges(At::Y(500.0), any_x)?;
    let footing = features.search_edges(At::Y(500.0), |x| x[0] <= 50.0)?;
    // let mut foot_ids = features.get_points_via_2d_edges(&footing);
    // let ids: Vec<_> = foot_ids.iter().map(|id| 1 + id).collect();
    // println!("Footing(IDs) = {:?}", ids);

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            // stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    /*
    const UY: [f64; 15] = [
        0.0,     //  0
        -0.0001,  //  1
        -0.00015, //  2
        -0.0002,  //  3
        -0.00025, //  4
        -0.00035, //  5
        -0.00045, //  6
        -0.00055, //  7
        -0.00065, //  8
        -0.00075, //  9
        -0.00080, // 10
        -0.00090, // 11
        -0.00110, // 12
        -0.00140, // 13
        -0.002,   // 14
    ];
    */

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&footing, Dof::Uy, -0.01);
    // println!("{}", essential);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_tol_rr(1e-6)
        .set_t_ini(0.0)
        .set_dt(|_| 1.0)
        .set_dt_out(|_| 1.0)
        .set_t_fin(1.0)
        // .set_constant_tangent(true)
        // .set_ignore_jacobian_symmetry(true)
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, NAME, None)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    /*
    // verify the results
    let tol_displacement = 1e-14;
    let tol_stress = 1e-8;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        "spo_754_footing_elast.json",
        tol_displacement,
        tol_stress,
        2,
    )?;
    assert!(all_good);
    */

    Ok(())
}
