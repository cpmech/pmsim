use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use russell_lab::*;

const NAME: &str = "spo_754_footing";
const DRAW_MESH_AND_EXIT: bool = false;
const VERBOSE_LEVEL: usize = 0;

const YOUNG: f64 = 1e7; // Young's modulus
const POISSON: f64 = 0.48; // Poisson's coefficient
const Z_INI: f64 = 848.7; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_754_footing() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", NAME))?;
    if DRAW_MESH_AND_EXIT {
        mesh.check_all()?;
        let mut fig = Figure::new();
        return fig
            .size(800.0, 800.0)
            .zoom_2d(15.0, 69.0, 448.0, 502.0, 0.5, 0.5, 0.5, 0.5)
            .range_2d(-10.0, 600.0, -10.0, 600.0)
            .draw(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME));
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let right = features.search_edges(At::X(500.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
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

    // displacement control
    const UY: [f64; 15] = [
        0.0,    // 0   t =  0.0
        -0.01,  // 1   t =  1.0
        -0.015, // 2   t =  2.0
        -0.02,  // 3   t =  3.0
        -0.025, // 4   t =  4.0
        -0.035, // 5   t =  5.0
        -0.045, // 6   t =  6.0
        -0.055, // 7   t =  7.0
        -0.065, // 8   t =  8.0
        -0.075, // 9   t =  9.0
        -0.08,  // 10  t = 10.0
        -0.09,  // 11  t = 11.0
        -0.11,  // 12  t = 12.0
        -0.14,  // 13  t = 13.0
        -0.2,   // 14  t = 14.0
    ];

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&right, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&footing, Dof::Uy, |t| UY[t as usize]);
    // println!("{}", essential);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_incremental(UY.len())
        // .set_constant_tangent(true)
        // .set_ignore_jacobian_symmetry(true)
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

    // verify the results
    let tol_displacement = 1e-10;
    let tol_stress = 5e-5;
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
    Ok(())
}
