use gemlab::mesh::{At, Cell, Extract, Mesh, Point, Region};
use gemlab::shapes::GeoKind;
use pmsim::base::{Config, Dof, FnBc, Nbc, ParamElement, ParamSolid, ParamStressStrain};
use pmsim::StrError;

fn main() -> Result<(), StrError> {
    // mesh
    //      y
    //      ^
    // 1.0  3------------2
    //      |`.      [1] |    [#] indicates id
    //      |  `.    (1) |    (#) indicates attribute_id
    //      |    `.      |
    //      |      `.    |
    //      | [0]    `.  |
    //      | (1)      `.|
    // 0.0  0------------1 -> x
    //     0.0          1.0
    #[rustfmt::skip]
    let mesh = Mesh {
        ndim: 2,
        points: vec![
            Point { id: 0, coords: vec![0.0, 0.0] },
            Point { id: 1, coords: vec![1.0, 0.0] },
            Point { id: 2, coords: vec![1.0, 1.0] },
            Point { id: 3, coords: vec![0.0, 1.0] },
        ],
        cells: vec![
            Cell { id: 0, attribute_id: 1, kind: GeoKind::Tri3, points: vec![0, 1, 3] },
            Cell { id: 1, attribute_id: 1, kind: GeoKind::Tri3, points: vec![2, 3, 1] },
        ],
    };

    // region and boundaries
    let region = Region::with(&mesh, Extract::Boundary)?;
    let bottom = region.find.edges(At::Y(0.0))?;
    let left = region.find.edges(At::X(0.0))?;
    let top = region.find.edges(At::Y(1.0))?;

    // configuration
    let mut config = Config::new(&region);
    let qn: FnBc = |_, _, _| -1.0;
    config
        .ebc_edges(&bottom, &[Dof::Uy], Config::zero)?
        .ebc_edges(&left, &[Dof::Ux], Config::zero)?
        .nbc_edges(&top, &[Nbc::Qn], qn)?;

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
         ndim = 2\n\
         npoint = 4\n\
         ncell = 2\n\
         nedge = 4\n\
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
         (1, Uy) @ t=0 → 0.0 @ t=1 → 0.0\n\
         (3, Ux) @ t=0 → 0.0 @ t=1 → 0.0\n\
         \n\
         Point boundary conditions\n\
         =========================\n\
         \n\
         Natural boundary conditions at edges\n\
         ====================================\n\
         ((2, 3), Qn) @ t=0 → -1.0 @ t=1 → -1.0\n\
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
