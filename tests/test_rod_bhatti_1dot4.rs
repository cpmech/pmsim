use pmsim::base::SampleMeshes;
use pmsim::fem::{sim_transient, BoundaryPointVec};
use pmsim::prelude::*;
use pmsim::StrError;
use russell_chk::vec_approx_eq;

// Bhatti's Example 1.4 example
#[test]
fn test_rod_bhatti_1dot4() -> Result<(), StrError> {
    //               (3)
    //               [2]
    //     2----------------------3
    //     |'.  (4)           _.-'
    //     |  '.[3]       _.-'
    //     |    '.    _.-'  (1)
    // (2) |      '1-'      [1]
    // [2] |      /
    //     |     /
    //     |    / (0)   The lines are ROD (Lin2) elements
    //     |   /  [1]
    //     |  /
    //     | /    (#) indicates cell id
    //     0'     [#] indicates attribute id
    let mesh = SampleMeshes::bhatti_example_1dot4_truss();

    // parameters, DOFs, and configuration
    #[rustfmt::skip]
    let data = Data::new(&mesh, [
        (1, Element::Rod(ParamRod { area: 4_000.0, young: 200_000.0, density: 1.0 })),
        (2, Element::Rod(ParamRod { area: 3_000.0, young: 200_000.0, density: 1.0 })),
        (3, Element::Rod(ParamRod { area: 2_000.0, young:  70_000.0, density: 1.0 })),
    ])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.at(&[0, 3], &[Dof::Ux, Dof::Uy], |_| 0.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.at(&[1], Pbc::Fy(|_| -150000.0));

    // boundary points
    let boundary_points = BoundaryPointVec::new(&data, &natural)?;

    // boundary elements
    let mut boundary_elements = BoundaryElementVec::new(&data, &config, &natural)?;

    // interior elements
    let mut interior_elements = InteriorElementVec::new(&data, &config)?;

    // simulation state
    let mut state = State::new(&data, &config, &essential)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).unwrap();

    // run simulation
    sim_transient(
        Some(&boundary_points),
        &mut boundary_elements,
        &mut interior_elements,
        &mut state,
        &mut lin_sys,
        &config,
    )?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00,  0.000000000000000e+00, // 0: Ux,Uy
        5.389536380057675e-01, -9.530613006371175e-01, // 1: Ux,Uy
        2.647036149579491e-01, -2.647036149579491e-01, // 2: Ux,Uy
        0.000000000000000e+00,  0.000000000000000e+00, // 3: Ux,Uy
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-15);
    Ok(())
}
