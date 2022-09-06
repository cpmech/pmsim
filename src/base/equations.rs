use super::{Dof, ElementInfoMap};
use crate::StrError;
use gemlab::mesh::{Mesh, PointId};
use std::collections::{HashMap, HashSet};
use std::fmt;

/// Holds equation numbers (DOF numbers)
///
/// # Examples
///
/// ## Given the mesh:
///
/// ```text
/// leq: local equation number       leq   point   geq
/// geq: global equation number       ↓        ↓    ↓
///                                   0 → Ux @ 0 →  0
///            {Ux → 6}               1 → Uy @ 0 →  1
///            {Uy → 7}               2 → Ux @ 1 →  3
///            {Pl → 8}               3 → Uy @ 1 →  4
///                2                  4 → Ux @ 2 →  6
///               / \                 5 → Uy @ 2 →  7
///   {Ux → 13}  /   \  {Ux → 11}     6 → Ux @ 3 →  9
///   {Uy → 14} 5     4 {Uy → 12}     7 → Uy @ 3 → 10
///            /       \              8 → Ux @ 4 → 11
/// {Ux → 0}  /         \  {Ux → 3}   9 → Uy @ 4 → 12
/// {Uy → 1} 0-----3-----1 {Uy → 4}  10 → Ux @ 5 → 13
/// {Pl → 2}   {Ux → 9}    {Pl → 5}  11 → Uy @ 5 → 14
///            {Uy → 10}             12 → Pl @ 0 →  2  <<< eq_first_pl
///                                  13 → Pl @ 1 →  5
///                                  14 → Pl @ 2 →  8
/// ```
///
/// ## Print the DOF numbers:
///
/// ```
/// use gemlab::mesh::Samples;
/// use gemlab::StrError;
/// use pmsim::base::{Attributes, Dof, Equations, Element, ElementInfoMap, SampleParams};
/// use std::collections::HashMap;
///
/// fn main() -> Result<(), StrError> {
///     let mesh = Samples::one_tri6();
///     let p1 = SampleParams::param_porous_sld_liq();
///     let att = Attributes::from([(1, Element::PorousSldLiq(p1))]);
///     let emap = ElementInfoMap::new(&mesh, &att)?;
///     let mut eqs = Equations::new(&mesh, &emap)?;
///     assert_eq!(
///         format!("{}", eqs),
/// r#"Points: DOFs and global equation numbers
/// ========================================
/// 0: [(Ux, 0), (Uy, 1), (Pl, 2)]
/// 1: [(Ux, 3), (Uy, 4), (Pl, 5)]
/// 2: [(Ux, 6), (Uy, 7), (Pl, 8)]
/// 3: [(Ux, 9), (Uy, 10)]
/// 4: [(Ux, 11), (Uy, 12)]
/// 5: [(Ux, 13), (Uy, 14)]
/// "#
///     );
///     Ok(())
/// }
/// ```
pub struct Equations {
    /// Holds all points DOFs and numbers
    ///
    /// **Notes:**
    ///
    /// 1. The array has a length equal to npoint
    /// 2. The inner maps have variable lengths according to the number of DOFs at the point
    pub all: Vec<HashMap<Dof, usize>>,

    /// Holds the total number of global equations
    ///
    /// **Note:** This is equal to the total number of DOFs
    pub n_equation: usize,
}

impl Equations {
    /// Allocates a new instance
    pub fn new(mesh: &Mesh, emap: &ElementInfoMap) -> Result<Self, StrError> {
        // auxiliary memoization data
        let npoint = mesh.points.len();
        let mut memo_point_dofs = vec![HashSet::new(); npoint];

        // find all element DOFs and local numbers and add (unique) DOF numbers to the point DOFs array
        for cell in &mesh.cells {
            let info = emap.get(cell)?;
            for m in 0..cell.points.len() {
                for (dof, _) in &info.dofs[m] {
                    memo_point_dofs[cell.points[m]].insert(*dof);
                }
            }
        }

        // compute all point DOF numbers
        let mut all = vec![HashMap::new(); npoint];
        let mut n_equation = 0; // equals the total number of DOFs
        for point_id in 0..npoint {
            let mut sorted_dofs: Vec<_> = memo_point_dofs[point_id].iter().collect();
            sorted_dofs.sort();
            for dof in sorted_dofs {
                all[point_id].insert(*dof, n_equation);
                n_equation += 1;
            }
        }

        // done
        Ok(Equations { all, n_equation })
    }

    /// Returns the (global) equation number of a (PointId,DOF) pair
    pub fn eq(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        if point_id >= self.all.len() {
            return Err("cannot find equation number because PointId is out-of-bounds");
        }
        let eq = self.all[point_id]
            .get(&dof)
            .ok_or("cannot find equation number corresponding to (PointId,DOF)")?;
        Ok(*eq)
    }
}

impl fmt::Display for Equations {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Points: DOFs and global equation numbers\n").unwrap();
        write!(f, "========================================\n").unwrap();
        for point_id in 0..self.all.len() {
            let mut dof_eqn: Vec<_> = self.all[point_id].iter().collect();
            dof_eqn.sort_by(|a, b| a.0.partial_cmp(b.0).unwrap());
            write!(f, "{:?}: {:?}\n", point_id, dof_eqn).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Equations;
    use crate::base::{Attributes, Dof, Element, ElementInfoMap, SampleParams};
    use gemlab::mesh::{PointId, Samples};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri6();
        let mut mesh_wrong = mesh.clone();
        mesh_wrong.cells[0].attribute_id = 100; // << never do this!
        let p1 = SampleParams::param_solid();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        assert_eq!(
            Equations::new(&mesh_wrong, &emap).err(),
            Some("cannot find (CellAttributeId, GeoKind) in ElementInfoMap")
        );
    }

    fn assert_equations(eqs: &Equations, p: PointId, correct: &[(Dof, usize)]) {
        let mut dofs: Vec<_> = eqs.all[p].iter().map(|(d, n)| (*d, *n)).collect();
        dofs.sort();
        assert_eq!(dofs, correct);
    }

    #[test]
    fn new_and_eq_works() {
        //                   {15}
        //     {21}          {16}
        //     {22}   {19}   {17}
        //     {23}   {20}   {18}  {13}
        //      8------7------6._  {14}
        //      |       [3](3)|  '-.5
        //      |  [0]        |     '-._
        // {24} 9  (1)      *10  [1]    '4 {11}
        // {25} |             |  (2)  .-'  {12}
        //      |       [2](3)|   _.3'
        //      0------1------2.-' {9}
        //     {0}    {3}    {5}   {10}
        //     {1}    {4}    {6}
        //     {2}           {7}
        //                   {8}  *10 => {26,27,28}
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_beam();
        let att = Attributes::from([
            (1, Element::PorousSldLiq(p1)),
            (2, Element::Solid(p2)),
            (3, Element::Beam(p3)),
        ]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();

        // check point dofs
        assert_equations(&eqs, 0, &[(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 2)]);
        assert_equations(&eqs, 1, &[(Dof::Ux, 3), (Dof::Uy, 4)]);
        assert_equations(&eqs, 2, &[(Dof::Ux, 5), (Dof::Uy, 6), (Dof::Rz, 7), (Dof::Pl, 8)]);
        assert_equations(&eqs, 3, &[(Dof::Ux, 9), (Dof::Uy, 10)]);
        assert_equations(&eqs, 4, &[(Dof::Ux, 11), (Dof::Uy, 12)]);
        assert_equations(&eqs, 5, &[(Dof::Ux, 13), (Dof::Uy, 14)]);
        assert_equations(&eqs, 6, &[(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Rz, 17), (Dof::Pl, 18)]);
        assert_equations(&eqs, 7, &[(Dof::Ux, 19), (Dof::Uy, 20)]);
        assert_equations(&eqs, 8, &[(Dof::Ux, 21), (Dof::Uy, 22), (Dof::Pl, 23)]);
        assert_equations(&eqs, 9, &[(Dof::Ux, 24), (Dof::Uy, 25)]);
        assert_equations(&eqs, 10, &[(Dof::Ux, 26), (Dof::Uy, 27), (Dof::Rz, 28)]);

        // check eq
        assert_eq!(eqs.eq(0, Dof::Ux).unwrap(), 0);
        assert_eq!(eqs.eq(6, Dof::Rz).unwrap(), 17);
        assert_eq!(eqs.eq(8, Dof::Pl).unwrap(), 23);
        assert_eq!(
            eqs.eq(111, Dof::Ux).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        assert_eq!(
            eqs.eq(0, Dof::T).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );
    }

    #[test]
    fn display_works() {
        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute_id
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_solid();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
        assert_eq!(
            format!("{}", eqs),
            "Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n"
        );

        // 3------------2------------5
        // |`.      [1] |            |    [#] indicates id
        // |  `.    (1) |            |    (#) indicates attribute_id
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let p = SampleParams::param_porous_liq();
        let att = Attributes::from([(1, Element::PorousLiq(p)), (2, Element::PorousLiq(p))]);
        let emap = ElementInfoMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
        assert_eq!(
            format!("{}", eqs),
            "Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Pl, 0)]\n\
             1: [(Pl, 1)]\n\
             2: [(Pl, 2)]\n\
             3: [(Pl, 3)]\n\
             4: [(Pl, 4)]\n\
             5: [(Pl, 5)]\n"
        );
    }
}
