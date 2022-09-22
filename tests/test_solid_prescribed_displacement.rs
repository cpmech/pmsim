use gemlab::mesh::Samples;
use pmsim::base::assemble_matrix;
use pmsim::fem::{ElementSolid, LocalEquations, PrescribedValues};
use pmsim::{prelude::*, StrError};
use russell_chk::vec_approx_eq;
use russell_lab::prelude::*;
use russell_sparse::prelude::*;

#[test]
fn test_solid_prescribed_displacement_direct() -> Result<(), StrError> {
    //  Uy PRESCRIBED          Uy PRESCRIBED
    //    {Ux → 6 ❓}          {Ux → 4 ❓}
    //    {Uy → 7 ✅}          {Uy → 5 ✅}
    //             3------------2
    //             |            |     PLANE-STRAIN
    //             |            |     E = 1
    //             |            |     nu = 0.25
    //             |            |
    // {Ux → 0 ✅} |            | {Ux → 2 ❓}
    // {Uy → 1 ✅} 0------------1 {Uy → 3 ✅}
    //           Ux,Uy          Uy
    //           FIXED         FIXED
    let mesh = Samples::one_qua4();

    // parameters
    const YOUNG: f64 = 1.0;
    const POISSON: f64 = 0.25;
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
    };

    // prescribed strain value (compression)
    const PRESSURE: f64 = -2.0;
    const STRAIN_Y: f64 = -PRESSURE * (POISSON - 1.0) * (POISSON + 1.0) / YOUNG;
    let (q, ee, nu) = (PRESSURE, YOUNG, POISSON);
    let strain = vec![-q * nu * (nu + 1.0) / ee, -q * (nu - 1.0) * (nu + 1.0) / ee, 0.0, 0.0];
    let stress = vec![0.0, q, q * nu, 0.0];

    // data, DOF numbers, and equations
    let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
    println!("\n{}", data.equations);
    let neq = data.equations.n_equation;
    assert_eq!(neq, 8);

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .at(&[0], Ebc::Ux(|_| 0.0))
        .at(&[0, 1], Ebc::Uy(|_| 0.0))
        .at(&[2, 3], Ebc::Uy(|_| STRAIN_Y));
    let values = PrescribedValues::new(&data, &essential)?;
    let prescribed = &values.flags;
    println!("                 0     1      2     3      4     5      6     7");
    println!("prescribed = {:?}", prescribed);

    // prescribed and unknown equations
    let mut eq_prescribed = values.equations.clone();
    let eq_unknown: Vec<_> = (0..neq).into_iter().filter(|i| !prescribed[*i]).collect();
    eq_prescribed.sort();
    println!("eq_prescribed = {:?}", eq_prescribed);
    println!("eq_unknown    = {:?}", eq_unknown);
    assert_eq!(eq_prescribed, &[0, 1, 3, 5, 7]);
    assert_eq!(eq_unknown, &[2, 4, 6]);

    // element and state
    let config = Config::new();
    let mut elem = ElementSolid::new(&data, &config, &mesh.cells[0], &p1).unwrap();
    let mut state = State::new(&data, &config)?;

    // update state with prescribed displacements
    let mut duu = Vector::new(neq); // cumulated ΔU
    values.apply(&mut duu, &mut state.uu, 1.0);
    println!("\nU = \n{:.6}", state.uu);
    println!("ΔU = \n{:.6}", duu);

    // local Jacobian matrix
    let neq_local = neq; // because of one element only
    let mut kk_local = Matrix::new(neq_local, neq_local);
    elem.calc_jacobian(&mut kk_local, &state)?;

    // full global Jacobian matrix
    let mut kk_global = SparseTriplet::new(neq, neq * neq)?;
    assemble_matrix(&mut kk_global, &kk_local, &elem.local_to_global, &vec![false; neq]);
    let mut kk_global_mat = Matrix::new(neq, neq);
    kk_global.to_matrix(&mut kk_global_mat)?;
    println!("\nK =\n{:.6}", kk_global_mat);

    // partitioned system
    let (n_prescribed, n_unknown) = (eq_prescribed.len(), eq_unknown.len());
    let mut kk11 = Matrix::new(n_unknown, n_unknown);
    let mut kk12 = Matrix::new(n_unknown, n_prescribed);
    for (i, ii) in eq_unknown.iter().enumerate() {
        for (j, jj) in eq_unknown.iter().enumerate() {
            kk11[i][j] = kk_global_mat[*ii][*jj];
        }
        for (j, jj) in eq_prescribed.iter().enumerate() {
            kk12[i][j] = kk_global_mat[*ii][*jj];
        }
    }
    println!("K11 = \n{:.6}", kk11);
    println!("K12 = \n{:.6}", kk12);
    let mut uu1 = Vector::new(n_unknown); // unknown displacements vector {U1}
    let mut uu2 = Vector::new(n_prescribed); // prescribed displacements vector {U2}
    let mut ee1 = Vector::new(n_unknown); // external forces vector {E1}
    for (i, ii) in eq_prescribed.iter().enumerate() {
        uu2[i] = state.uu[*ii];
    }
    println!("U2 = \n{:.6}", uu2);

    // fix vector of external forces
    mat_vec_mul(&mut ee1, -1.0, &kk12, &uu2)?;
    println!("E1 = \n{:.6}", ee1);

    // solve linear system
    vec_copy(&mut uu1, &ee1)?;
    solve_lin_sys(&mut uu1, &mut kk11)?;
    println!("U1 = \n{}", uu1);

    // fix vector of displacements and increments
    for (i, ii) in eq_unknown.iter().enumerate() {
        state.uu[*ii] = uu1[i];
        duu[*ii] = uu1[i];
    }
    println!("\nU = \n{:.6}", state.uu);
    println!("ΔU = \n{:.6}", duu);

    // update stresses
    println!("\nstrain_y = {:?}", STRAIN_Y);
    println!("strain = {:?}", strain);
    println!("stress = {:?}", stress);
    elem.update_secondary_values(&state, &duu)?;
    for p in 0..elem.ips.len() {
        println!("σ = {:?}", elem.stresses[p].sigma.vec.as_data());
        vec_approx_eq(elem.stresses[p].sigma.vec.as_data(), &stress, 1e-15);
    }

    // compute external forces
    let mut kk21 = Matrix::new(n_prescribed, n_unknown);
    let mut kk22 = Matrix::new(n_prescribed, n_prescribed);
    for (i, ii) in eq_prescribed.iter().enumerate() {
        for (j, jj) in eq_unknown.iter().enumerate() {
            kk21[i][j] = kk_global_mat[*ii][*jj];
        }
        for (j, jj) in eq_prescribed.iter().enumerate() {
            kk22[i][j] = kk_global_mat[*ii][*jj];
        }
    }
    println!("K21 = \n{:.6}", kk21);
    println!("K22 = \n{:.6}", kk22);
    let mut ee2 = Vector::new(n_prescribed); // external forces vector {E2} = [K21]{U1} + [K22]{U2}
    let mut ee2_tmp = Vector::new(n_prescribed);
    mat_vec_mul(&mut ee2, 1.0, &kk21, &uu1)?;
    mat_vec_mul(&mut ee2_tmp, 1.0, &kk22, &uu2)?;
    vec_update(&mut ee2, 1.0, &ee2_tmp)?;
    println!("F_ext(2) = \n{}", ee2);

    // compute residual (actually, internal forces)
    let mut rr_local = Vector::new(neq_local);
    elem.calc_residual(&mut rr_local, &state)?;
    println!("\nF_int = \n{}", rr_local);
    Ok(())
}
