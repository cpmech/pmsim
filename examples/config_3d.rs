use gemlab::mesh::{At, Cell, Extract, Mesh, Point, Region};
use gemlab::shapes::GeoKind;
use pmsim::base::{Config, Dof, FnBc, Nbc, ParamElement, ParamSolid, ParamStressStrain};
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    // mesh
    //       4--------------7  1.0
    //      /.             /|
    //     / .            / |    [#] indicates id
    //    /  .           /  |    (#) indicates attribute_id
    //   /   .          /   |
    //  5--------------6    |          z
    //  |    .         |    |          ↑
    //  |    0---------|----3  0.0     o → y
    //  |   /  [0]     |   /          ↙
    //  |  /   (1)     |  /          x
    //  | /            | /
    //  |/             |/
    //  1--------------2   1.0
    // 0.0            1.0
    #[rustfmt::skip]
    let mesh = Mesh {
        ndim: 3,
        points: vec![
            Point { id: 0, coords: vec![0.0, 0.0, 0.0] },
            Point { id: 1, coords: vec![1.0, 0.0, 0.0] },
            Point { id: 2, coords: vec![1.0, 1.0, 0.0] },
            Point { id: 3, coords: vec![0.0, 1.0, 0.0] },
            Point { id: 4, coords: vec![0.0, 0.0, 1.0] },
            Point { id: 5, coords: vec![1.0, 0.0, 1.0] },
            Point { id: 6, coords: vec![1.0, 1.0, 1.0] },
            Point { id: 7, coords: vec![0.0, 1.0, 1.0] },
        ],
        cells: vec![
            Cell { id: 0, attribute_id: 1, kind: GeoKind::Hex8, points: vec![0,1,2,3, 4,5,6,7] },
        ],
    };

    // region and boundaries
    let region = Region::with(&mesh, Extract::Boundary)?;
    let x_zero = region.find.faces(At::X(0.0))?;
    let y_zero = region.find.faces(At::Y(0.0))?;
    let z_zero = region.find.faces(At::Z(0.0))?;
    let top = region.find.faces(At::Z(1.0))?;

    // configuration
    let mut config = Config::new(&region);
    let qn: FnBc = |_, _, _| -1.0;
    config
        .ebc_faces(&x_zero, &[Dof::Ux], Config::zero)?
        .ebc_faces(&y_zero, &[Dof::Uy], Config::zero)?
        .ebc_faces(&z_zero, &[Dof::Uz], Config::zero)?
        .nbc_faces(&top, &[Nbc::Qn], qn)?;

    // parameters and elements
    let solid = ParamSolid {
        density: 2.7, // Mg/m²
        stress_strain: ParamStressStrain::LinearElastic {
            young: 10_000.0, // kPa
            poisson: 0.2,    // [-]
        },
        n_integ_point: None,
    };
    config.elements(1, ParamElement::Solid(solid))?;

    assert_eq!(
        format!("{}", config),
        "Mesh data\n\
         =========\n\
         ndim = 3\n\
         npoint = 8\n\
         ncell = 1\n\
         nedge = 12\n\
         nface = 6\n\
         \n\
         Other configuration data\n\
         ========================\n\
         gravity = 0.0\n\
         thickness = 1.0\n\
         plane_stress = false\n\
         total_stress = false\n\
         initialization = Zero\n\
         \n\
         Essential boundary conditions\n\
         =============================\n\
         (0, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (0, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (0, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (1, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (1, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (2, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (3, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (3, Uz) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (4, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (4, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (5, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (7, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
         \n\
         Point boundary conditions\n\
         =========================\n\
         \n\
         Natural boundary conditions at edges\n\
         ====================================\n\
         \n\
         Natural boundary conditions at faces\n\
         ====================================\n\
         ((4, 5, 6, 7), Qn) @ t=0 → -1.0 @ t=1 → -1.0\n\
         \n\
         Parameters for fluids\n\
         =====================\n\
         None\n\
         \n\
         Parameters for Elements\n\
         =======================\n\
         1 → Solid(ParamSolid { density: 2.7, stress_strain: LinearElastic { young: 10000.0, poisson: 0.2 }, n_integ_point: None })\n"
        );
    Ok(())
}
