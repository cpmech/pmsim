#![allow(unused)]

use gemlab::mesh::{EdgeKey, Extract, Features, PointId};
use gemlab::shapes::{GeoKind, Scratchpad};
use pmsim::base::{
    assemble_matrix, assemble_vector, BcEssential, BcNatural, Dof, DofNumbers, Element, FnBc, Nbc, ParamSolid,
    ParamStressStrain, SampleMeshes,
};
use pmsim::StrError;
use russell_chk::assert_vec_approx_eq;
use russell_lab::{Matrix, Vector};
use russell_sparse::{ConfigSolver, Solver};
use russell_sparse::{SparseTriplet, Symmetry};
use std::collections::{HashMap, HashSet};

/*
fn linear_solution(dn: &DofNumbers, bc: &Conditions, elements: &[Solid]) -> Result<(), StrError> {
    // allocate coefficient matrix and rhs vector
    let neq = dn.n_equation;
    let nnz = dn.nnz_sup;
    let mut kk_global = SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap();
    let mut rr_global = Vector::new(neq);

    let mut pads: Vec<Scratchpad> = Vec::new();
    pads.push(Scratchpad::new(2, GeoKind::Tri3)?);

    // assembly
    for i in 0..elements.len() {
        let e = &elements[i];
        assemble_matrix(&mut kk_global, &e.kk_local, &dn.local_to_global[i], &dn.prescribed);
        assemble_vector(&mut rr_global, &e.r_local, &dn.local_to_global[i], &dn.prescribed);
    }
    for eq in 0..neq {
        if dn.prescribed[eq] {
            kk_global.put(eq, eq, 1.0)?;
        }
    }

    // solve linear system
    let config = ConfigSolver::new();
    let (_, uu) = Solver::compute(config, &kk_global, &rr_global)?;
    println!("{}", uu);
    Ok(())
}
*/

// Solution of Bhatti's Example 1.6
#[test]
fn test_rod_bhatti_1dot6() -> Result<(), StrError> {
    // 2.0  fixed 1'-,_load                connectivity:
    //            |     '-,_      load      eid : vertices
    // 1.5 - - -  |        ,'3-,__            0 :  0, 2, 3
    //            |  1   ,'  |    '-,_        1 :  3, 1, 0
    // 1.0 - - -  |    ,'    |  3   ,-'5      2 :  2, 4, 5
    //            |  ,'  0   |   ,-'   |      3 :  5, 3, 2
    //            |,'        |,-'   2  |
    // 0.0  fixed 0----------2---------4   constraints:
    //           0.0        2.0       4.0   fixed on x and y
    let mesh = SampleMeshes::bhatti_example_1dot6_bracket();
    let features = Features::new(&mesh, Extract::Boundary);

    // DOF numbers
    let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid)]))?;

    // parameters
    let param = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 10_000.0,
            poisson: 0.2,
        },
        two_dim: true,
        plane_stress: true,
        thickness: 0.25,
        n_integ_point: None,
    };

    let mut ebc = BcEssential::new();
    let zero = |_| 0.0;
    ebc.points(&[0, 1], &[Dof::Ux, Dof::Uy], zero);

    let mut nbc = BcNatural::new();
    nbc.edge_keys(&features, &[(1, 3), (3, 5)], Nbc::Qn, |_| -20.0)?;
    // nbcs.edges(
    //     // &[&feat.edges.get(&(1, 3)).unwrap(), &feat.edges.get(&(3, 5)).unwrap()],
    //     feat.get_edges([(1, 3), (3, 5)]).unwrap(),
    //     Nbc::Qn,
    //     |_| -20.0,
    // );

    // let mut bc = Conditions::new();
    // bc.set_essential_at_points(&region, &HashSet::from([0, 1]), &[Dof::Ux, Dof::Uy], zero)?;
    // bc.set_natural_at_edges(&region, &HashSet::from([(1, 3), (3, 5)]), Nbc::Qn, |_, _, _| -20.0)?;

    // elements
    // let mut elements = Vec::new();
    // for cell in &mesh.cells {
    //     elements.push(Solid::new(&mesh, cell.id, &param, &dn, plane_stress, thickness)?);
    // }

    // // solution
    // linear_solution(&dn, &bc, &elements);
    Ok(())
}
