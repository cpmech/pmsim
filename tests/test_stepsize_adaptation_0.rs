use gemlab::prelude::*;
use pmsim::prelude::*;
use pmsim::StrError;

// This test catches the failure due to large δu
//
// See test_spo_751_pres_cylin for details on the problem.

const NAME_MESH: &str = "spo_751_pres_cylin";

const A: f64 = 100.0; // inner radius

const P_MAX_RES: f64 = 0.20; // maximum pressure applied before unloading completely to zero
const PP: [f64; 2] = [
    P_MAX_RES, // stage = 0
    0.0,       // stage = 1
];

const YOUNG: f64 = 210.0; // Young's modulus
const POISSON: f64 = 0.3; // Poisson's coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_stepsize_adaptation_0() -> Result<(), StrError> {
    // mesh
    let kind = GeoKind::Qua4;
    let mesh = Mesh::read(&format!("data/spo/{}_{}.msh", NAME_MESH, kind.to_string())).unwrap();

    // features
    let features = Features::new(&mesh, false);
    let bottom = features.search_edges(At::Y(0.0), any_x)?;
    let left = features.search_edges(At::X(0.0), any_x)?;
    let inner_circle = features.search_edges(At::Circle(0.0, 0.0, A), any_x)?;

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

    // configuration
    let mut config = Config::new(&mesh);
    config.set_steady(PP.len()).set_substepping(false);

    // natural boundary conditions and configuration
    let mut natural = Natural::new();
    natural.edges_fn(&inner_circle, Nbc::Qn, |stage, _| -PP[stage]);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new();

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // check if the solver indeed failed
    assert_eq!(solver.has_failed(), true);
    Ok(())
}
