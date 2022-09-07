use pmsim::base::SampleMeshes;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_chk::vec_approx_eq;

// Bhatti's Example 1.4 on page 25
#[test]
fn test_rod_bhatti_1dot4() -> Result<(), StrError> {
    // mesh and boundary features
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
    let mut config = Config::new();
    config.control.n_max_time_steps = 2;

    // essential boundary conditions
    let mut essential = Essential::new();
    let zero = |_| 0.0;
    essential.at(&[0, 3], Ebc::Ux(zero)).at(&[0, 3], Ebc::Uy(zero));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.at(&[1], Pbc::Fy(|_| -150000.0));

    // point loads
    let concentrated_loads = ConcentratedLoads::new(&data, &natural)?;

    // prescribed values
    let prescribed_values = PrescribedValues::new(&data, &essential)?;

    // boundary elements
    let mut boundary_elements = BoundaryElements::new(&data, &config, &natural)?;

    // interior elements
    let mut interior_elements = InteriorElements::new(&data, &config)?;

    // simulation state
    let mut state = State::new(&data, &config)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &prescribed_values, &interior_elements, &boundary_elements)?;

    // run simulation
    simulation(
        Some(&concentrated_loads),
        &prescribed_values,
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
