use super::{Attributes, ElementDofs};
use crate::StrError;
use gemlab::mesh::{Cell, CellAttribute, Mesh};
use gemlab::shapes::GeoKind;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;

/// Maps (CellAttribute, GeoKind) to ElementDofs
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ElementDofsMap {
    all: HashMap<(CellAttribute, GeoKind), ElementDofs>,
}

impl ElementDofsMap {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, att_map: &Attributes) -> Result<Self, StrError> {
        let mut all = HashMap::new();
        for cell in &mesh.cells {
            let element = att_map.get(cell.attribute)?;
            all.insert(
                (cell.attribute, cell.kind),
                ElementDofs::new(mesh.ndim, *element, cell.kind)?,
            );
        }
        Ok(ElementDofsMap { all })
    }

    /// Returns the ElementDofs corresponding to Cell
    pub fn get(&self, cell: &Cell) -> Result<&ElementDofs, StrError> {
        self.all
            .get(&(cell.attribute, cell.kind))
            .ok_or("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
    }
}

impl fmt::Display for ElementDofsMap {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Elements: DOFs and local equation numbers\n").unwrap();
        write!(f, "=========================================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort_by(|a, b| a.0.cmp(&b.0));
        for key in keys {
            let info = self.all.get(key).unwrap();
            let (id, kind) = key;
            write!(f, "{} → {:?}\n", id, kind).unwrap();
            write!(f, "{}", info).unwrap();
            write!(f, "-----------------------------------------\n").unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDofsMap;
    use crate::base::{Attributes, Elem, ParamPorousLiq};
    use crate::base::{ParamBeam, ParamPorousSldLiq, ParamRod, ParamSolid};
    use gemlab::mesh::Samples;

    #[test]
    fn new_map_handles_errors() {
        let mesh = Samples::one_tri6();
        let p2 = ParamSolid::sample_linear_elastic();
        let amap = Attributes::from([(2, Elem::Solid(p2))]);
        assert_eq!(
            ElementDofsMap::new(&mesh, &amap).err(),
            Some("cannot find CellAttribute in Attributes map")
        );
        let p1 = ParamRod::sample();
        let amap = Attributes::from([(1, Elem::Rod(p1))]);
        assert_eq!(
            ElementDofsMap::new(&mesh, &amap).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    #[test]
    fn new_map_and_get_work() {
        let mesh = Samples::three_tri3();
        let mut mesh_wrong = mesh.clone();
        let p1 = ParamSolid::sample_linear_elastic();
        let amap = Attributes::from([(1, Elem::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        assert_eq!(emap.get(&mesh.cells[0]).unwrap().n_equation, 6);
        mesh_wrong.cells[0].attribute = 100; // never do this
        assert_eq!(
            emap.get(&mesh_wrong.cells[0]).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
    }

    #[test]
    fn new_map_display_works() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let amap = Attributes::from([(1, Elem::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        assert_eq!(
            format!("{}", emap),
            "Elements: DOFs and local equation numbers\n\
             =========================================\n\
             1 → Tri3\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n"
        );

        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let amap = Attributes::from([(1, Elem::PorousLiq(p)), (2, Elem::PorousLiq(p))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        assert_eq!(
            format!("{}", emap),
            "Elements: DOFs and local equation numbers\n\
             =========================================\n\
             1 → Tri3\n\
             0: [(Pl, 0)]\n\
             1: [(Pl, 1)]\n\
             2: [(Pl, 2)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n\
             2 → Qua4\n\
             0: [(Pl, 0)]\n\
             1: [(Pl, 1)]\n\
             2: [(Pl, 2)]\n\
             3: [(Pl, 3)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n"
        );

        // 8------7------6._
        // |       [3](3)|  '-.5
        // |  [0]        |     '-._
        // 9  (1)       10  [1]    '4
        // |             |  (2)  .-'
        // |       [2](3)|   _.3'
        // 0------1------2.-'
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let amap = Attributes::from([(1, Elem::PorousSldLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Beam(p3))]);
        let emap = ElementDofsMap::new(&mesh, &amap).unwrap();
        assert_eq!(
            format!("{}", emap),
            "Elements: DOFs and local equation numbers\n\
             =========================================\n\
             1 → Qua8\n\
             0: [(Ux, 0), (Uy, 1), (Pl, 16)]\n\
             1: [(Ux, 2), (Uy, 3), (Pl, 17)]\n\
             2: [(Ux, 4), (Uy, 5), (Pl, 18)]\n\
             3: [(Ux, 6), (Uy, 7), (Pl, 19)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             5: [(Ux, 10), (Uy, 11)]\n\
             6: [(Ux, 12), (Uy, 13)]\n\
             7: [(Ux, 14), (Uy, 15)]\n\
             (Pl @ Some(16), Pg @ None, T @ None)\n\
             -----------------------------------------\n\
             2 → Tri6\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             5: [(Ux, 10), (Uy, 11)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n\
             3 → Lin2\n\
             0: [(Ux, 0), (Uy, 1), (Rz, 2)]\n\
             1: [(Ux, 3), (Uy, 4), (Rz, 5)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n"
        );
    }
}
