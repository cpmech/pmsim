use super::{Dof, ElementDofsMap};
use crate::StrError;
use gemlab::mesh::{Mesh, PointId};
use serde::{Deserialize, Serialize};
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
/// use pmsim::base::{Attributes, Dof, Equations, ElementDofsMap, Etype};
/// use pmsim::base::ParamPorousSldLiq;
/// use std::collections::HashMap;
///
/// fn main() -> Result<(), StrError> {
///     let mesh = Samples::one_tri6();
///     let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
///     let att = Attributes::from([(1, Etype::PorousSldLiq(p1))]);
///     let emap = ElementDofsMap::new(&mesh, &att)?;
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
#[derive(Clone, Debug, Serialize, Deserialize)]
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
    pub fn new(mesh: &Mesh, emap: &ElementDofsMap) -> Result<Self, StrError> {
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
    use crate::base::{Attributes, Dof, ElementDofsMap, Etype, SampleMeshes};
    use crate::base::{ParamBeam, ParamPorousLiq, ParamPorousSldLiq, ParamSolid};
    use gemlab::mesh::{PointId, Samples};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri6();
        let mut mesh_wrong = mesh.clone();
        mesh_wrong.cells[0].attribute = 100; // << never do this!
        let p1 = ParamSolid::sample_linear_elastic();
        let att = Attributes::from([(1, Etype::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        assert_eq!(
            Equations::new(&mesh_wrong, &emap).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
    }

    fn assert_equations(eqs: &Equations, p: PointId, correct: &[(Dof, usize)]) {
        let mut dofs: Vec<_> = eqs.all[p].iter().map(|(d, n)| (*d, *n)).collect();
        dofs.sort();
        assert_eq!(dofs, correct);
    }

    #[test]
    fn new_and_eq_works() {
        //                     {Ux→15}
        //    {Ux→21}          {Uy→16}
        //    {Uy→22}  {Ux→19} {Rz→17}
        //    {Pl→23}  {Uy→20} {Pl→18} {Ux→13}
        //         8------7------6._   {Uy→14}
        //         |       [3](3)|  '-.5
        //         |  [0]        |     '-._
        // {Ux→24} 9  (1)      *10  [1]    '4 {Ux→11}
        // {Uy→25} |             |  (2)  .-'  {Uy→12}
        //         |       [2](3)|   _.3'
        //         0------1------2.-'  {Ux→9}
        //     {Ux→0}  {Ux→3}  {Ux→5}  {Uy→10}
        //     {Uy→1}  {Uy→4}  {Uy→6}
        //     {Pl→2}          {Rz→7}
        //                     {Pl→8}
        //  *10 => {Ux→26, Uy→27, Rz→28}
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let att = Attributes::from([
            (1, Etype::PorousSldLiq(p1)),
            (2, Etype::Solid(p2)),
            (3, Etype::Beam(p3)),
        ]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
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
        let att = Attributes::from([(1, Etype::Solid(p1))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
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
        // |  `.    (1) |            |    (#) indicates attribute
        // |    `.      |     [2]    |
        // |      `.    |     (2)    |
        // | [0]    `.  |            |
        // | (1)      `.|            |
        // 0------------1------------4
        let mesh = Samples::two_tri3_one_qua4();
        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let att = Attributes::from([(1, Etype::PorousLiq(p)), (2, Etype::PorousLiq(p))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
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

    #[test]
    fn coupled_dofs_work() {
        //                   {Ux→38}
        //        {Ux→24}    {Uy→39}    {Ux→27}
        //        {Uy→25} 8----14-----9 {Uy→28}
        //        {Pl→26} |           | {Pl→29}
        //                |           |
        // {Ux→52,Uy→53} 21    26    22 {Ux→54,Uy→55}   26: {Ux→62,Uy→63}
        //                |           |
        //        {Ux→18} |  {Ux→36}  | {Ux→21}
        //        {Uy→19} 6----13-----7 {Uy→22}
        //        {Pl→20} |  {Uy→37}  | {Pl→23}
        //                |           |
        // {Ux→48,Uy→49} 19    25    20 {Ux→50,Uy→51}   25: {Ux→60,Uy→61}
        //                |           |
        //        {Ux→12} |  {Ux→34}  | {Ux→15}
        //        {Uy→13} 4----12-----5 {Uy→16}
        //        {Pl→14} |  {Uy→35}  | {Pl→17}
        //                |           |
        // {Ux→44,Uy→45} 17    24    18 {Ux→46,Uy→47}   24: {Ux→58,Uy→59}
        //                |           |
        //         {Ux→6} |  {Ux→32}  | {Ux→9}
        //         {Uy→7} 2----11-----3 {Uy→10}
        //         {Pl→8} |  {Uy→33}  | {Pl→11}
        //                |           |
        // {Ux→40,Uy→41} 15    23    16 {Ux→42,Uy→43}   23: {Ux→56,Uy→57}
        //                |           |
        //         {Ux→0} |           | {Ux→3}
        //         {Uy→1} 0----10-----1 {Uy→4}
        //         {Pl→2}    {Ux→30}    {Pl→5}
        //                   {Uy→31}
        let mesh = SampleMeshes::column_two_layers_qua9();
        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let att = Attributes::from([(1, Etype::PorousSldLiq(p)), (2, Etype::PorousSldLiq(p))]);
        let emap = ElementDofsMap::new(&mesh, &att).unwrap();
        let eqs = Equations::new(&mesh, &emap).unwrap();
        assert_eq!(
            format!("{}", eqs),
            "Points: DOFs and global equation numbers\n\
             ========================================\n\
             0: [(Ux, 0), (Uy, 1), (Pl, 2)]\n\
             1: [(Ux, 3), (Uy, 4), (Pl, 5)]\n\
             2: [(Ux, 6), (Uy, 7), (Pl, 8)]\n\
             3: [(Ux, 9), (Uy, 10), (Pl, 11)]\n\
             4: [(Ux, 12), (Uy, 13), (Pl, 14)]\n\
             5: [(Ux, 15), (Uy, 16), (Pl, 17)]\n\
             6: [(Ux, 18), (Uy, 19), (Pl, 20)]\n\
             7: [(Ux, 21), (Uy, 22), (Pl, 23)]\n\
             8: [(Ux, 24), (Uy, 25), (Pl, 26)]\n\
             9: [(Ux, 27), (Uy, 28), (Pl, 29)]\n\
             10: [(Ux, 30), (Uy, 31)]\n\
             11: [(Ux, 32), (Uy, 33)]\n\
             12: [(Ux, 34), (Uy, 35)]\n\
             13: [(Ux, 36), (Uy, 37)]\n\
             14: [(Ux, 38), (Uy, 39)]\n\
             15: [(Ux, 40), (Uy, 41)]\n\
             16: [(Ux, 42), (Uy, 43)]\n\
             17: [(Ux, 44), (Uy, 45)]\n\
             18: [(Ux, 46), (Uy, 47)]\n\
             19: [(Ux, 48), (Uy, 49)]\n\
             20: [(Ux, 50), (Uy, 51)]\n\
             21: [(Ux, 52), (Uy, 53)]\n\
             22: [(Ux, 54), (Uy, 55)]\n\
             23: [(Ux, 56), (Uy, 57)]\n\
             24: [(Ux, 58), (Uy, 59)]\n\
             25: [(Ux, 60), (Uy, 61)]\n\
             26: [(Ux, 62), (Uy, 63)]\n"
        );
    }
}
