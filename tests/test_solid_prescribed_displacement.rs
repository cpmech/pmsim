use gemlab::mesh::Samples;
use pmsim::base::{assemble_matrix, assemble_vector};
use pmsim::fem::{BcPrescribedArray, ElementSolid, ElementTrait};
use pmsim::prelude::*;
use russell_lab::*;
use russell_sparse::prelude::*;

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

const YOUNG: f64 = 1.0;
const POISSON: f64 = 0.25;
const PRESSURE: f64 = -2.0;
const STRAIN_Y: f64 = -PRESSURE * (POISSON - 1.0) * (POISSON + 1.0) / YOUNG;

#[test]
fn test_solid_prescribed_displacement_direct_approach() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
        ngauss: None,
    };

    // prescribed strain value (compression)
    let (q, ee, nu) = (PRESSURE, YOUNG, POISSON);
    let strain = vec![-q * nu * (nu + 1.0) / ee, -q * (nu - 1.0) * (nu + 1.0) / ee, 0.0, 0.0];
    let stress = vec![0.0, q, q * nu, 0.0];

    // data, DOF numbers, and equations
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
    let neq = base.dofs.size();
    assert_eq!(neq, 8);

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .points(&[0], Dof::Ux, 0.0)
        .points(&[0, 1], Dof::Uy, 0.0)
        .points(&[2, 3], Dof::Uy, STRAIN_Y);
    let values = BcPrescribedArray::new(&base, &essential)?;
    let prescribed = &values.flags;

    // prescribed and unknown equations
    let mut eq_prescribed = values.equations.clone();
    let eq_unknown: Vec<_> = (0..neq).into_iter().filter(|i| !prescribed[*i]).collect();
    eq_prescribed.sort();
    assert_eq!(eq_prescribed, &[0, 1, 3, 5, 7]);
    assert_eq!(eq_unknown, &[2, 4, 6]);

    // element and state
    let config = Config::new(&mesh);
    let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // pub fn apply(&self, duu: &mut Vector, uu: &mut Vector, time: f64) {
    // self.all.iter().for_each(|e| e.set_value(duu, uu, time));
    // }

    // update state with prescribed displacements
    for p in 0..values.equations.len() {
        let value = values.all[p].value(1.0);
        let eq = values.equations[p];
        state.duu[eq] = value - state.uu[eq];
        state.uu[eq] = value;
    }
    println!("\nU = \n{:.6}", state.uu);
    println!("ΔU = \n{:.6}", state.duu);

    // global = local Jacobian matrix (kk_global = kk_local because there is one element only)
    let mut kk = Matrix::new(neq, neq);
    elem.calc_jacobian(&mut kk, &state)?;

    // partitioned system
    let (n_prescribed, n_unknown) = (eq_prescribed.len(), eq_unknown.len());
    let mut kk11 = Matrix::new(n_unknown, n_unknown);
    let mut kk12 = Matrix::new(n_unknown, n_prescribed);
    for (i, ii) in eq_unknown.iter().enumerate() {
        for (j, jj) in eq_unknown.iter().enumerate() {
            kk11.set(i, j, kk.get(*ii, *jj));
        }
        for (j, jj) in eq_prescribed.iter().enumerate() {
            kk12.set(i, j, kk.get(*ii, *jj));
        }
    }
    let mut uu1 = Vector::new(n_unknown); // unknown displacements vector U1
    let mut uu2 = Vector::new(n_prescribed); // prescribed displacements vector U2
    let mut ee1 = Vector::new(n_unknown); // external forces vector E1
    for (i, ii) in eq_prescribed.iter().enumerate() {
        uu2[i] = state.uu[*ii];
    }

    // fix vector of external forces: E1 = -K12·U2
    mat_vec_mul(&mut ee1, -1.0, &kk12, &uu2)?;
    println!("E1 = \n{:.6}", ee1);

    // solve linear system
    vec_copy(&mut uu1, &ee1)?; // RHS = U1 ← E1 (just for the solver)
    solve_lin_sys(&mut uu1, &mut kk11)?;
    // println!("U1 = \n{}", uu1);

    // fix vector of displacements and increments
    for (i, ii) in eq_unknown.iter().enumerate() {
        state.uu[*ii] = uu1[i];
        state.duu[*ii] = uu1[i];
    }
    println!("\nU = \n{:.6}", state.uu);
    println!("ΔU = \n{:.6}", state.duu);

    // update stresses
    println!("\nstrain = {:?}", strain);
    println!("stress = {:?}", stress);
    elem.update_secondary_values(&mut state)?;
    for p in 0..elem.gauss.npoint() {
        println!("σ = {:?}", state.gauss[0].stress(p)?.vector().as_data());
        vec_approx_eq(&state.gauss[0].stress(p)?.vector(), &stress, 1e-15);
    }

    // compute external forces
    let mut kk21 = Matrix::new(n_prescribed, n_unknown);
    let mut kk22 = Matrix::new(n_prescribed, n_prescribed);
    for (i, ii) in eq_prescribed.iter().enumerate() {
        for (j, jj) in eq_unknown.iter().enumerate() {
            kk21.set(i, j, kk.get(*ii, *jj));
        }
        for (j, jj) in eq_prescribed.iter().enumerate() {
            kk22.set(i, j, kk.get(*ii, *jj));
        }
    }
    let mut ee2 = Vector::new(n_prescribed); // external forces vector E2 = K21·U1 + K22·U2
    let mut ee2_tmp = Vector::new(n_prescribed);
    mat_vec_mul(&mut ee2, 1.0, &kk21, &uu1)?;
    mat_vec_mul(&mut ee2_tmp, 1.0, &kk22, &uu2)?;
    vec_update(&mut ee2, 1.0, &ee2_tmp)?;
    println!("\nexternal F2 = \n{}", ee2);

    // compute (local=global) residual (actually, internal forces)
    let mut rr = Vector::new(neq);
    elem.calc_residual(&mut rr, &state)?;
    let mut rr2 = Vector::new(n_prescribed);
    for (i, ii) in eq_prescribed.iter().enumerate() {
        rr2[i] = rr[*ii];
    }
    println!("internal F2 = \n{}", rr2);
    vec_approx_eq(&ee2, &rr2, 1e-15);
    Ok(())
}

#[test]
fn test_solid_prescribed_displacement_residual_approach() -> Result<(), StrError> {
    // mesh
    let mesh = Samples::one_qua4();

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
        ngauss: None,
    };

    // prescribed strain value (compression)
    let (q, ee, nu) = (PRESSURE, YOUNG, POISSON);
    let strain = vec![-q * nu * (nu + 1.0) / ee, -q * (nu - 1.0) * (nu + 1.0) / ee, 0.0, 0.0];
    let stress = vec![0.0, q, q * nu, 0.0];

    // data, DOF numbers, and equations
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
    let neq = base.dofs.size();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .points(&[0], Dof::Ux, 0.0)
        .points(&[0, 1], Dof::Uy, 0.0)
        .points(&[2, 3], Dof::Uy, STRAIN_Y);
    let values = BcPrescribedArray::new(&base, &essential)?;
    let prescribed = &values.flags;

    // prescribed and unknown equations
    let eq_prescribed = values.equations.clone();
    let eq_unknown: Vec<_> = (0..neq).into_iter().filter(|i| !prescribed[*i]).collect();
    let n_unknown = eq_unknown.len();

    // element and state
    let config = Config::new(&mesh);
    let mut elem = ElementSolid::new(&mesh, &base, &config, &p1, 0).unwrap();
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // update state with prescribed displacements
    for p in 0..values.equations.len() {
        let value = values.all[p].value(1.0);
        let eq = values.equations[p];
        state.duu[eq] = value - state.uu[eq];
        state.uu[eq] = value;
    }
    println!("\nU = \n{:.6}", state.uu);
    println!("ΔU = \n{:.6}", state.duu);

    // update secondary variables (corresponds to E1 = -K12·U2)
    println!("\nstrain = {:?}", strain);
    println!("stress = {:?}", stress);
    elem.update_secondary_values(&mut state)?;
    for p in 0..elem.gauss.npoint() {
        println!("σ = {:?}", state.gauss[0].stress(p)?.vector().as_data());
    }

    // compute residual (actually, internal forces)
    let mut rr_local = Vector::new(neq);
    elem.calc_residual(&mut rr_local, &state)?;
    println!("local R = \n{}", rr_local);
    let mut ee1 = Vector::new(n_unknown); // external forces vector E1
    for (i, ii) in eq_unknown.iter().enumerate() {
        ee1[i] = -rr_local[*ii];
    }
    assert_eq!(
        format!("{:.6}", ee1),
        "┌           ┐\n\
         │  0.375000 │\n\
         │  0.375000 │\n\
         │ -0.375000 │\n\
         └           ┘"
    );

    // assemble modified residual vector
    let mut rr_global = Vector::new(neq);
    assemble_vector(&mut rr_global, &rr_local, &elem.local_to_global, &prescribed);
    println!("modified global R = \n{}", rr_global);

    // local Jacobian matrix (== global Jacobian matrix)
    let mut kk_local = Matrix::new(neq, neq);
    elem.calc_jacobian(&mut kk_local, &state)?;

    // global Jacobian matrix
    let tol = Some(1e-14);
    let mut kk_global = SparseMatrix::new_coo(neq, neq, neq * neq, Sym::No)?;
    assemble_matrix(
        kk_global.get_coo_mut()?,
        &kk_local,
        &elem.local_to_global,
        &prescribed,
        tol,
    )
    .unwrap();
    for eq in &values.equations {
        kk_global.put(*eq, *eq, 1.0)?;
    }

    // solve linear system
    let mut mdu = Vector::new(neq);
    LinSolver::compute(Genie::Umfpack, &mut mdu, &mut kk_global, &rr_global, None)?;

    // update U and ΔU
    for i in 0..neq {
        state.uu[i] -= mdu[i];
        state.duu[i] -= mdu[i];
    }
    let duu_expected = vec![0.0, 0.0, 0.625, 0.0, 0.625, -1.875, 0.0, -1.875];
    vec_approx_eq(&state.duu, &duu_expected, 1e-15);

    // must reset ΔU corresponding to prescribed values because
    // it has been used to update the stress already!
    for eq in &eq_prescribed {
        state.duu[*eq] = 0.0;
    }

    // update secondary variables
    println!("\nstrain = {:?}", strain);
    println!("stress = {:?}", stress);
    elem.update_secondary_values(&mut state)?;
    // σ = [0.0, -2.0000000000000004, -0.5, 4.4408920985006264e-17]
    for p in 0..elem.gauss.npoint() {
        println!("σ = {:?}", state.gauss[0].stress(p)?.vector().as_data());
        vec_approx_eq(&state.gauss[0].stress(p)?.vector(), &stress, 1e-15);
    }
    Ok(())
}
