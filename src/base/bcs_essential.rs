use super::{Dof, DofNumbers, FnBc};
use crate::StrError;
use gemlab::mesh::{Feature, PointId};
use std::collections::HashMap;
use std::fmt;

/// Holds essential boundary conditions
pub struct BcsEssential {
    pub all: HashMap<(PointId, Dof), FnBc>,
}

impl BcsEssential {
    /// Allocates a new instance
    pub fn new() -> Self {
        BcsEssential { all: HashMap::new() }
    }

    /// Sets essential boundary condition at points
    pub fn at(&mut self, points: &[PointId], dofs: &[Dof], f: FnBc) -> &mut Self {
        for point_id in points {
            for dof in dofs {
                self.all.insert((*point_id, *dof), f);
            }
        }
        self
    }

    /// Sets essential boundary condition on edges or faces
    pub fn on(&mut self, features: &[&Feature], dofs: &[Dof], f: FnBc) -> &mut Self {
        for edge in features {
            for point_id in &edge.points {
                for dof in dofs {
                    self.all.insert((*point_id, *dof), f);
                }
            }
        }
        self
    }

    /// Generates two arrays to handle prescribed DOFs
    ///
    /// # Output
    ///
    /// * `prescribed` -- Is an `n_equation` array of `bool` indicating which DOFs (equations) are prescribed.
    ///   The length of `prescribed` is equal to the total number of DOFs (total number of equations).
    /// * `p_equations` -- Is a "smaller" array with only the DOFs numbers of the prescribed equations.
    pub fn prescribed(&self, dn: &DofNumbers) -> Result<(Vec<bool>, Vec<usize>), StrError> {
        let mut prescribed = vec![false; dn.n_equation];
        let mut p_equations = Vec::new();
        for (point_id, dof) in self.all.keys() {
            match dn.point_dofs[*point_id].get(dof) {
                Some(eq) => {
                    prescribed[*eq] = true;
                    p_equations.push(*eq)
                }
                None => return Err("EBC dof is not present in point_dofs array"),
            }
        }
        Ok((prescribed, p_equations))
    }
}

impl fmt::Display for BcsEssential {
    /// Prints a formatted summary of Boundary Conditions
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Essential boundary conditions\n").unwrap();
        write!(f, "=============================\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort();
        for key in keys {
            let fbc = self.all.get(key).unwrap();
            let f0 = fbc(0.0);
            let f1 = fbc(1.0);
            write!(f, "{:?} : {:?} @ t=0 → {:?} @ t=1 → {:?}\n", key.0, key.1, f0, f1).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcsEssential;
    use crate::base::{Attributes, Dof, DofNumbers, Element, ElementDofsMap, SampleParams};
    use gemlab::mesh::{Feature, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn bcs_essential_works() {
        let mut bcs_essential = BcsEssential::new();
        let edges = &[&Feature {
            kind: GeoKind::Lin2,
            points: vec![1, 2],
        }];
        let faces = &[&Feature {
            kind: GeoKind::Tri3,
            points: vec![3, 4, 5],
        }];
        bcs_essential
            .at(&[0], &[Dof::Ux, Dof::Uy], |_| 0.0)
            .on(edges, &[Dof::Pl], |t| t)
            .on(faces, &[Dof::T], |t| t / 2.0);
        assert_eq!(
            format!("{}", bcs_essential),
            "Essential boundary conditions\n\
             =============================\n\
             0 : Ux @ t=0 → 0.0 @ t=1 → 0.0\n\
             0 : Uy @ t=0 → 0.0 @ t=1 → 0.0\n\
             1 : Pl @ t=0 → 0.0 @ t=1 → 1.0\n\
             2 : Pl @ t=0 → 0.0 @ t=1 → 1.0\n\
             3 : T @ t=0 → 0.0 @ t=1 → 0.5\n\
             4 : T @ t=0 → 0.0 @ t=1 → 0.5\n\
             5 : T @ t=0 → 0.0 @ t=1 → 0.5\n"
        );
    }

    #[test]
    fn prescribed_works() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute_id
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_porous_liq();
        let att = Attributes::from([(1, Element::PorousLiq(p1))]);
        let edm = ElementDofsMap::new(&mesh, &att).unwrap();
        let dn = DofNumbers::new(&mesh, &edm).unwrap();
        let mut bcs_essential = BcsEssential::new();
        let zero = |_| 0.0;
        assert_eq!(zero(1.0), 0.0);
        bcs_essential.at(&[0, 4], &[Dof::Pl], zero);
        let (prescribed, mut p_equations) = bcs_essential.prescribed(&dn).unwrap();
        assert_eq!(prescribed, &[true, false, false, false, true]);
        p_equations.sort();
        assert_eq!(p_equations, &[0, 4]);

        bcs_essential.at(&[3], &[Dof::T], zero);
        assert_eq!(
            bcs_essential.prescribed(&dn).err(),
            Some("EBC dof is not present in point_dofs array")
        );

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
        let p1 = SampleParams::param_solid();
        let att = Attributes::from([(1, Element::Solid(p1))]);
        let edm = ElementDofsMap::new(&mesh, &att).unwrap();
        let dn = DofNumbers::new(&mesh, &edm).unwrap();
        let mut bcs_essential = BcsEssential::new();
        bcs_essential.at(&[0], &[Dof::Ux, Dof::Uy], zero);
        bcs_essential.at(&[1, 2], &[Dof::Uy], zero);
        let (prescribed, mut p_equations) = bcs_essential.prescribed(&dn).unwrap();
        assert_eq!(
            prescribed,
            //   0     1      2     3      4     5      6      7      8      9
            &[true, true, false, true, false, true, false, false, false, false]
        );
        p_equations.sort();
        assert_eq!(p_equations, &[0, 1, 3, 5]);
    }
}
