#![allow(unused)]

use gemlab::mesh::{At, Find};
use pmsim::base::{Config, Dof, Element, Essential, Natural, Nbc, ParamDiffusion, SampleMeshes};
use pmsim::fem::{BoundaryElementVec, Data, InteriorElementVec, LinearSystem, State};
use pmsim::StrError;
use russell_chk::vec_approx_eq;
use russell_lab::{copy_vector, Matrix, Vector};

#[test]
fn test_bhatti_6dot22_heat() -> Result<(), StrError> {
    // mesh, parameters, DOFs, and configuration
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

    // boundary features
    let find = Find::new(&mesh, None); // boundary only
    let bottom = find.edges(At::Y(0.0))?;
    let edges_flux = find.edges(At::X(0.0))?;
    let edges_conv = vec![
        find.edges(At::Y(0.03))?.as_slice(),  // top-horizontal
        find.edges(At::X(0.03))?.as_slice(),  // middle-vertical
        find.edges(At::Y(0.015))?.as_slice(), // middle-horizontal
    ]
    .concat();
    println!("flux: {:?}", edges_flux);
    println!("conv: {:?}", edges_conv);

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

    // boundary elements
    let mut boundary_elements = BoundaryElementVec::new(&data, &config, &natural)?;

    // interior elements
    let mut interior_elements = InteriorElementVec::new(&data, &config)?;

    // simulation state
    let mut state = State::new(&data, &config, &essential)?;

    // check residual of first element
    state.uu.fill(0.0);
    // with state = 0, the residual is equal to -b (negative of integral of source term)
    let neg_b = Vector::from(&[250.0, 312.5, 312.5, 250.0, -1125., -1000.0, -1125., -1250.0]);
    let res: Result<(), _> = interior_elements
        .all
        .iter_mut()
        .map(|e| {
            e.calc_residual(&state)?;
            e.calc_jacobian(&state)
        })
        .collect();
    if let Some(err) = res.err() {
        return Err(err);
    }
    let elem0 = &interior_elements.all[0];
    let elem1 = &interior_elements.all[1];
    vec_approx_eq(elem0.residual.as_data(), neg_b.as_data(), 1e-12);

    // check Jacobian of first element (independent of state)
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
    vec_approx_eq(elem0.jacobian.as_data(), bhatti_kk0.as_data(), 1e-13);

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
    vec_approx_eq(elem1.jacobian.as_data(), bhatti_kk1.as_data(), 1e-12);

    // fix state.uu (must do this before calculating residuals)
    let (prescribed, p_equations) = data.prescribed(&essential)?;
    for eq in &p_equations {
        state.uu[*eq] = 110.0;
    }

    /*
    // assemble system
    let mut lin_sys = LinearSystem::new(&data);
    let rr = &mut lin_sys.residual;
    let kk = &mut lin_sys.jacobian;
    elements.iter_mut().for_each(|e| {
        e.calc_residual(&state)?;
        e.calc_jacobian(&state)?;
        assemble_vector(rr, &e.residual, &e.local_to_global, &prescribed);
        assemble_matrix(kk, &e.jacobian, &e.local_to_global, &prescribed);
    });
    nbcs.all.iter_mut().for_each(|e| {
        e.calc_residual(&state);
        e.calc_jacobian(&state);
        // println!("{}", e.residual);
        assemble_vector(rr, &e.residual, &e.local_to_global, &prescribed);
        match &e.jacobian {
            Some(jj) => {
                // println!("{}", jj);
                assemble_matrix(kk, jj, &e.local_to_global, &prescribed);
            }
            None => (),
        }
    });
    for i in 0..neq {
        if prescribed[i] {
            kk.put(i, i, 1.0);
        }
    }

    println!("rr =\n{}", rr);
    // let mut kk_mat = Matrix::new(neq, neq);
    // kk.to_matrix(&mut kk_mat);
    // println!("kk =\n{:.4}", kk_mat);

    let bhatti_rr = &[
        204.5, -1147.0, 304.25, -1011.0, 616.75, -1022.0, 307.0, -1125.0, 0.0, 0.0, 0.0, -1410.0, -2250.0,
    ];
    // vec_approx_eq(rr.as_data(), bhatti_rr, 1e-12);

    let mut mdu = Vector::new(neq);
    lin_sys.solver.initialize(&kk)?;
    lin_sys.solver.factorize()?;
    lin_sys.solver.solve(&mut mdu, &rr);

    let uu = &state.primary_unknowns;
    // println!("uu =\n{}", uu);
    let uu_new = &mut Vector::new(neq);
    add_vectors(uu_new, 1.0, uu, -1.0, &mdu);
    println!("uu_new =\n{}", uu_new);

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
    copy_vector(&mut state.primary_unknowns, &uu_new);

    rr.fill(0.0);
    elements.iter_mut().for_each(|e| {
        e.calc_residual(&state)?;
        assemble_vector(rr, &e.residual, &e.local_to_global, &prescribed);
    });
    nbcs.all.iter_mut().for_each(|e| {
        e.calc_residual(&state);
        assemble_vector(rr, &e.residual, &e.local_to_global, &prescribed);
    });
    println!("rr_new =\n{}", rr);
    */
    Ok(())
}
