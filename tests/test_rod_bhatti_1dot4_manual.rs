use pmsim::base::{assemble_matrix, assemble_vector, AttrElement, DataMaps, Element, ParamRod, SampleMeshes};
use pmsim::element::Rod;
use pmsim::StrError;
use russell_chk::assert_vec_approx_eq;
use russell_lab::{Matrix, Vector};
use russell_sparse::{ConfigSolver, Solver};
use russell_sparse::{SparseTriplet, Symmetry};
use std::collections::HashMap;

// Manual implementation of Bhatti's Example 1.4 example
#[test]
fn test_rod_bhatti_1dot4_manual() -> Result<(), StrError> {
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

    // parameters
    #[rustfmt::skip]
    let attr_param = HashMap::from([
        (1, ParamRod { area: 4_000.0, young: 200_000.0, density: 1.0 }),
        (2, ParamRod { area: 3_000.0, young: 200_000.0, density: 1.0 }),
        (3, ParamRod { area: 2_000.0, young:  70_000.0, density: 1.0 }),
    ]);

    // elements
    let mut elements = Vec::new();
    let mut attr_elem = AttrElement::new();
    for cell in &mesh.cells {
        elements.push(Rod::new(&mesh, cell.id, attr_param.get(&cell.attribute_id).unwrap()));
        attr_elem.insert(cell.attribute_id, Element::Rod);
    }

    // datamaps
    let dm = DataMaps::new(&mesh, attr_elem)?;
    let (neq, nnz) = dm.neq_nnz();

    // prescribed equations
    let prescribed = vec![false; neq];

    // assembly
    let mut kk = SparseTriplet::new(neq, neq, nnz, Symmetry::No).unwrap();
    for i in 0..elements.len() {
        let rod = &elements[i];
        assemble_matrix(&mut kk, &rod.ke, i, &dm.local_to_global, &prescribed);
    }
    let mut kk_mat = Matrix::new(neq, neq);
    kk.to_matrix(&mut kk_mat).unwrap();

    // check full matrix
    #[rustfmt::skip]
    let kk_correct = [
          3.2600217813448358e+04,  7.6067174898046171e+04, -3.2600217813448358e+04, -7.6067174898046171e+04,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00,  0.0000000000000000e+00, // 0
          7.6067174898046171e+04,  2.9749007476210775e+05, -7.6067174898046171e+04, -1.7749007476210775e+05,  0.0000000000000000e+00, -1.2000000000000000e+05,  0.0000000000000000e+00,  0.0000000000000000e+00, // 1
         -3.2600217813448358e+04, -7.6067174898046171e+04,  2.4308860903092835e+05,  1.1913603334072011e+05, -3.2998316455372231e+04,  3.2998316455372231e+04, -1.7749007476210775e+05, -7.6067174898046171e+04, // 2
         -7.6067174898046171e+04, -1.7749007476210775e+05,  1.1913603334072011e+05,  2.4308860903092835e+05,  3.2998316455372231e+04, -3.2998316455372231e+04, -7.6067174898046171e+04, -3.2600217813448358e+04, // 3
          0.0000000000000000e+00,  0.0000000000000000e+00, -3.2998316455372231e+04,  3.2998316455372231e+04,  1.5299831645537223e+05, -3.2998316455372231e+04, -1.2000000000000000e+05,  0.0000000000000000e+00, // 4
          0.0000000000000000e+00, -1.2000000000000000e+05,  3.2998316455372231e+04, -3.2998316455372231e+04, -3.2998316455372231e+04,  1.5299831645537223e+05,  0.0000000000000000e+00,  0.0000000000000000e+00, // 5
          0.0000000000000000e+00,  0.0000000000000000e+00, -1.7749007476210775e+05, -7.6067174898046171e+04, -1.2000000000000000e+05,  0.0000000000000000e+00,  2.9749007476210775e+05,  7.6067174898046171e+04, // 6
          0.0000000000000000e+00,  0.0000000000000000e+00, -7.6067174898046171e+04, -3.2600217813448358e+04,  0.0000000000000000e+00,  0.0000000000000000e+00,  7.6067174898046171e+04,  3.2600217813448358e+04, // 7
    ];
    assert_vec_approx_eq!(kk_mat.as_data(), kk_correct, 1e-10);

    // prescribed equations
    let mut prescribed = vec![false; neq];
    prescribed[dm.point_equations[0][0]] = true; // Ux
    prescribed[dm.point_equations[0][1]] = true; // Uy
    prescribed[dm.point_equations[3][0]] = true; // Ux
    prescribed[dm.point_equations[3][1]] = true; // Uy

    // assembly
    kk.reset();
    let mut ff = Vector::new(neq);
    for i in 0..elements.len() {
        let rod = &elements[i];
        assemble_matrix(&mut kk, &rod.ke, i, &dm.local_to_global, &prescribed);
        assemble_vector(&mut ff, &rod.fe, i, &dm.local_to_global, &prescribed);
    }
    for i in 0..neq {
        if prescribed[i] {
            kk.put(i, i, 1.0)?;
        }
    }

    // point loads
    ff[dm.point_equations[1][1]] = -150000.0;

    // solve linear system
    let config = ConfigSolver::new();
    let (_, uu) = Solver::compute(config, &kk, &ff)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = [
        0.000000000000000e+00,  0.000000000000000e+00, // 0: Ux,Uy
        5.389536380057675e-01, -9.530613006371175e-01, // 1: Ux,Uy
        2.647036149579491e-01, -2.647036149579491e-01, // 2: Ux,Uy
        0.000000000000000e+00,  0.000000000000000e+00, // 3: Ux,Uy
    ];
    assert_vec_approx_eq!(uu.as_data(), uu_correct, 1e-15);
    Ok(())
}
