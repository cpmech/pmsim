use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, fem::sim_transient, prelude::*, StrError};
use russell_chk::vec_approx_eq;
use russell_lab::{add_vectors, copy_vector, mat_approx_eq, vector_norm, Matrix, NormVec, Vector};

#[test]
fn test_bhatti_6dot22_heat() -> Result<(), StrError> {
    // mesh and boundary features
    //       0.0    0.015    0.03
    // 0.03   0-------1-------2
    //        |               |
    //        |               3
    //        |               |
    // 0.015 11            _.'4-------5-------6 0.015
    //        |        _.-'                   |
    //        |    _.-12                      7 0.0075
    //        |_.-'                           |
    // 0.0   10---------------9---------------8 0.0
    //       0.0             0.03            0.06
    let mesh = SampleMeshes::bhatti_example_6dot22_heat();
    let find = Find::new(&mesh, None); // boundary only
    let bottom = find.edges(At::Y(0.0))?;
    let edges_flux = find.edges(At::X(0.0))?;
    let edges_conv = vec![
        find.edges(At::Y(0.03))?.as_slice(),  // top-horizontal
        find.edges(At::X(0.03))?.as_slice(),  // middle-vertical
        find.edges(At::Y(0.015))?.as_slice(), // middle-horizontal
    ]
    .concat();
    let points_flux: Vec<_> = edges_flux.iter().map(|f| &f.points).collect();
    let points_conv: Vec<_> = edges_conv.iter().map(|f| &f.points).collect();
    println!("flux: {:?}", points_flux);
    println!("conv: {:?}", points_conv);
    assert_eq!(points_flux[0], &[10, 0, 11]);
    assert_eq!(points_conv[0], &[0, 2, 1]);
    assert_eq!(points_conv[1], &[2, 4, 3]);
    assert_eq!(points_conv[2], &[4, 6, 5]);

    // parameters, DOFs, and configuration
    let (kx, ky) = (45.0, 45.0);
    let source = 5e6;
    let p1 = ParamDiffusion {
        rho: 0.0,
        kx,
        ky,
        kz: 0.0,
        source: Some(source),
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&bottom, &[Dof::T], |_| 110.0);
    println!("\n{}", essential);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .on(&edges_flux, Nbc::Qt(|_| 8000.0))
        .on(&edges_conv, Nbc::Cv(55.0, |_| 20.0));
    println!("{}", natural);

    // interior elements
    let mut interior_elements = InteriorElementVec::new(&data, &config)?;

    // boundary elements
    let mut boundary_elements = BoundaryElementVec::new(&data, &config, &natural)?;

    // simulation state
    let mut state = State::new(&data, &config, &essential)?;

    // check residual of first element
    state.uu.fill(0.0);
    // with state = 0, the residual is equal to -b (negative of integral of source term)
    let neg_b = Vector::from(&[250.0, 312.5, 312.5, 250.0, -1125., -1000.0, -1125., -1250.0]);
    interior_elements.calc_residuals(&state)?;
    vec_approx_eq(interior_elements.all[0].residual.as_data(), neg_b.as_data(), 1e-12);

    // check Jacobian of first element (independent of state)
    interior_elements.calc_jacobians(&state)?;
    #[rustfmt::skip]
    let bhatti_kk0 = Matrix::from(&[
        [38.515873015873005  , 23.194444444444443  , 21.46825396825396   , 22.02777777777777   , -20.317460317460306 , -30.91269841269842  , -14.682539682539685 , -39.293650793650784],
        [23.19444444444444   , 84.08730158730152   , 33.6111111111111    , 30.456349206349195  , -26.984126984126974 , -82.69841269841265  , -28.015873015873005 , -33.650793650793645],
        [21.46825396825396   , 33.6111111111111    , 58.313492063492056  , 23.19444444444444   , -19.36507936507935  , -68.17460317460318  , -20.63492063492062  , -28.412698412698408],
        [22.02777777777777   , 30.456349206349195  , 23.19444444444444   , 57.14682539682538   , -16.507936507936492 , -36.150793650793645 , -33.49206349206349  , -46.67460317460315 ],
        [-20.317460317460302 , -26.984126984126966 , -19.365079365079353 , -16.507936507936492 , 95.07936507936506   , 16.349206349206334  , -5.0793650793650915 , -23.174603174603188],
        [-30.91269841269841  , -82.69841269841265  , -68.17460317460318  , -36.15079365079365  , 16.34920634920633   , 155.87301587301585  , 3.650793650793659   , 42.06349206349206  ],
        [-14.682539682539685 , -28.015873015873005 , -20.63492063492062  , -33.49206349206349  , -5.079365079365089  , 3.650793650793652   , 95.07936507936506   , 3.174603174603193  ],
        [-39.293650793650784 , -33.65079365079365  , -28.412698412698408 , -46.67460317460315  , -23.174603174603188 , 42.06349206349206   , 3.1746031746031935  , 125.96825396825392 ],
    ]);
    mat_approx_eq(&interior_elements.all[0].jacobian, &bhatti_kk0, 1e-13);

    // check Jacobian of second element (independent of state)
    #[rustfmt::skip]
    let bhatti_kk1 = Matrix::from(&[
        [43.05158730158727   , 25.313492063492063  , 26.646825396825385  , 49.623015873015845  , 16.634920634920615  , -48.015873015873005 , -21.26984126984126  , -91.98412698412692 ], 
        [25.313492063492063  , 117.57539682539678  , 49.62301587301586   , 62.599206349206334  , -12.888888888888902 , -144.68253968253964 , -42.22222222222222  , -55.3174603174603  ], 
        [26.646825396825385  , 49.623015873015845  , 70.45634920634916   , 52.00396825396825   , -11.269841269841278 , -103.96825396825392 , -27.460317460317462 , -56.031746031745996], 
        [49.623015873015845  , 62.599206349206334  , 52.00396825396825   , 173.5515873015872   , -32.222222222222214 , -90.63492063492062  , -85.55555555555557  , -129.36507936507928], 
        [16.634920634920615  , -12.888888888888909 , -11.269841269841281 , -32.222222222222214 , 156.25396825396825  , 12.698412698412731  , -36.50793650793649  , -92.69841269841268 ], 
        [-48.015873015873005 , -144.6825396825396  , -103.96825396825393 , -90.63492063492062  , 12.698412698412724  , 258.8888888888888   , 14.603174603174601  , 101.11111111111106 ], 
        [-21.269841269841265 , -42.22222222222221  , -27.460317460317462 , -85.55555555555557  , -36.50793650793649  , 14.60317460317459   , 133.01587301587304  , 65.39682539682536  ], 
        [-91.98412698412692  , -55.3174603174603   , -56.031746031745996 , -129.36507936507928 , -92.6984126984127   , 101.11111111111107  , 65.39682539682536   , 258.8888888888888  ], 
    ]);
    mat_approx_eq(&interior_elements.all[1].jacobian, &bhatti_kk1, 1e-12);

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).unwrap();

    // fix state.uu (must do this before calculating residuals)
    for eq in &lin_sys.p_equations {
        state.uu[*eq] = 110.0;
    }

    // compute residuals in parallel
    interior_elements.calc_residuals_parallel(&state)?;
    boundary_elements.calc_residuals_parallel(&state)?;

    // assemble residuals
    let rr = &mut lin_sys.residual;
    interior_elements.assemble_residuals(rr, &lin_sys.prescribed);
    boundary_elements.assemble_residuals(rr, &lin_sys.prescribed);
    println!("rr =\n{}", rr);
    let bhatti_rr = &[
        2627.5555555555547,
        -2762.079365079365,
        2665.757936507936,
        -4411.396825396825,
        11968.138888888885,
        -12021.999999999995,
        7456.9999999999945,
        -20924.99999999999,
        0.0,
        0.0,
        0.0,
        -5732.301587301586,
        -30884.92063492062,
    ];
    vec_approx_eq(rr.as_data(), bhatti_rr, 1e-10);
    let norm_rr = vector_norm(rr, NormVec::Max);
    println!("norm_rr = {:?}", norm_rr);

    // compute jacobians in parallel
    interior_elements.calc_jacobians_parallel(&state)?;
    boundary_elements.calc_jacobians_parallel(&state)?;

    // assemble jacobians matrices
    let kk = &mut lin_sys.jacobian;
    interior_elements.assemble_jacobians(kk, &lin_sys.prescribed);
    boundary_elements.assemble_jacobians(kk, &lin_sys.prescribed);
    let mut kk_mat = Matrix::new(lin_sys.n_equation, lin_sys.n_equation);
    kk.to_matrix(&mut kk_mat)?;
    // println!("kk =\n{:.4}", kk_mat);

    // check global Jacobian matrix
    #[rustfmt::skip]
    let bhatti_kk = Matrix::from(&[
        [57.36682539682538   , -33.38206349206349  , 23.13944444444444   , -36.150793650793645 , 30.456349206349195  , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , -46.67460317460315  , -16.507936507936492],
        [-33.38206349206349  , 95.95936507936506   , -20.52492063492062  , 3.650793650793652   , -28.015873015873005 , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 3.174603174603193   , -5.079365079365089 ],
        [23.13944444444444   , -20.52492063492062  , 58.643492063492054  , -68.11960317460317  , 33.5836111111111    , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , -28.412698412698408 , -19.36507936507935 ],
        [-36.15079365079365  , 3.650793650793659   , -68.11960317460317  , 156.31301587301584  , -82.64341269841265  , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 42.06349206349206   , 16.34920634920633  ],
        [30.456349206349195  , -28.015873015873005 , 33.5836111111111    , -82.64341269841265  , 257.9688888888887   , -85.44555555555557  , 51.94896825396825   , -90.63492063492062  , 0.0 , 0.0, 0.0 , -33.650793650793645 , -156.34920634920627],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , -85.44555555555557  , 133.89587301587304  , -27.350317460317463 , 14.60317460317459   , 0.0 , 0.0, 0.0 , 0.0                 , 65.39682539682536  ],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , 51.94896825396825   , -27.350317460317463 , 70.67634920634916   , -103.96825396825392 , 0.0 , 0.0, 0.0 , 0.0                 , -56.031746031745996],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , -90.63492063492062  , 14.603174603174601  , -103.96825396825393 , 258.8888888888888   , 0.0 , 0.0, 0.0 , 0.0                 , 101.11111111111106 ],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 0.0                 , 0.0                ],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 0.0                 , 0.0                ],
        [0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 0.0                 , 0.0                ],
        [-46.67460317460315  , 3.1746031746031935  , -28.412698412698408 , 42.06349206349206   , -33.65079365079365  , 0.0                 , 0.0                 , 0.0                 , 0.0 , 0.0, 0.0 , 125.96825396825392  , -23.174603174603188],
        [-16.507936507936492 , -5.0793650793650915 , -19.365079365079353 , 16.349206349206334  , -156.34920634920624 , 65.39682539682536   , -56.031746031745996 , 101.11111111111107  , 0.0 , 0.0, 0.0 , -23.174603174603188 , 353.96825396825386 ],
    ]);
    mat_approx_eq(&kk_mat, &bhatti_kk, 1e-12);

    // augment global Jacobian matrix
    for eq in &lin_sys.p_equations {
        kk.put(*eq, *eq, 1.0)?;
    }

    // solve linear system
    let mdu = &mut lin_sys.mdu;
    lin_sys.solver.initialize(&kk)?;
    lin_sys.solver.factorize()?;
    lin_sys.solver.solve(mdu, &rr)?;

    // update U vector
    let mut uu_new = Vector::new(lin_sys.n_equation);
    add_vectors(&mut uu_new, 1.0, &state.uu, -1.0, &mdu)?;
    println!("uu_new =\n{}", uu_new);

    // check U vector
    let tt_bhatti = Vector::from(&[
        156.440502466202,
        150.75605418729847,
        149.19646294563637,
        144.2245542836661,
        133.8432701060946,
        124.00195294431063,
        121.74635727622194,
        119.14813150652589,
        110.0,
        110.0,
        110.0,
        144.67542222443012,
        129.13200798820264,
    ]);
    vec_approx_eq(uu_new.as_data(), tt_bhatti.as_data(), 1e-12);

    // set state with new U vector and check the residuals
    copy_vector(&mut state.uu, &uu_new)?;
    rr.fill(0.0);
    interior_elements.calc_residuals(&state)?;
    boundary_elements.calc_residuals(&state)?;
    interior_elements.assemble_residuals(rr, &lin_sys.prescribed);
    boundary_elements.assemble_residuals(rr, &lin_sys.prescribed);
    println!("rr_new =\n{:?}", rr);
    let norm_rr = vector_norm(rr, NormVec::Max);
    println!("norm_rr = {:?}", norm_rr);
    assert!(norm_rr < 1e-10);
    Ok(())
}

#[test]
fn test_bhatti_6dot22_heat_sim() -> Result<(), StrError> {
    // mesh and boundary features
    let mesh = SampleMeshes::bhatti_example_6dot22_heat();
    let find = Find::new(&mesh, None); // boundary only
    let bottom = find.edges(At::Y(0.0))?;
    let edges_flux = find.edges(At::X(0.0))?;
    let edges_conv = vec![
        find.edges(At::Y(0.03))?.as_slice(),  // top-horizontal
        find.edges(At::X(0.03))?.as_slice(),  // middle-vertical
        find.edges(At::Y(0.015))?.as_slice(), // middle-horizontal
    ]
    .concat();

    // parameters, DOFs, and configuration
    let (kx, ky) = (45.0, 45.0);
    let source = 5e6;
    let p1 = ParamDiffusion {
        rho: 1.0,
        kx,
        ky,
        kz: 0.0,
        source: Some(source),
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let mut config = Config::new();
    config.transient = true;
    config.control.n_max_time_steps = 2;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&bottom, &[Dof::T], |_| 110.0);

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .on(&edges_flux, Nbc::Qt(|_| 8000.0))
        .on(&edges_conv, Nbc::Cv(55.0, |_| 20.0));

    // interior elements
    let mut interior_elements = InteriorElementVec::new(&data, &config)?;

    // boundary elements
    let mut boundary_elements = BoundaryElementVec::new(&data, &config, &natural)?;

    // simulation state
    let mut state = State::new(&data, &config, &essential)?;

    // linear system
    let mut lin_sys = LinearSystem::new(&data, &essential, &interior_elements, &boundary_elements).unwrap();

    // run simulation
    sim_transient(
        &mut interior_elements,
        &mut boundary_elements,
        &mut state,
        &mut lin_sys,
        &config,
    )
}
