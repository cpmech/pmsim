use crate::base::{Attributes, Element, ElementInfoMap, Equations, Essential};
use crate::StrError;
use gemlab::mesh::{Cell, CellAttributeId, Mesh};

/// Data holds the basic data for a FEM simulation
pub struct Data<'a> {
    /// Holds an access to the Mesh
    pub mesh: &'a Mesh,

    /// Holds all attributes
    pub attributes: Attributes,

    /// Holds the element information such as local DOFs and equation numbers
    pub information: ElementInfoMap,

    /// Holds all DOF numbers
    pub equations: Equations,
}

impl<'a> Data<'a> {
    /// Allocate new instance
    pub fn new<const N: usize>(mesh: &'a Mesh, arr: [(CellAttributeId, Element); N]) -> Result<Self, StrError> {
        let attributes = Attributes::from(arr);
        let information = ElementInfoMap::new(&mesh, &attributes)?;
        let equations = Equations::new(&mesh, &information).unwrap(); // cannot fail
        Ok(Data {
            mesh,
            attributes,
            information,
            equations,
        })
    }

    /// Returns the number of local equations
    pub fn n_local_eq(&self, cell: &Cell) -> Result<usize, StrError> {
        let info = self.information.get(cell)?;
        Ok(info.n_equation)
    }

    /// Generates two arrays to handle prescribed DOFs
    ///
    /// # Output
    ///
    /// * `prescribed` -- Is an `n_equation` array of `bool` indicating which DOFs (equations) are prescribed.
    ///   The length of `prescribed` is equal to the total number of DOFs (total number of equations).
    /// * `p_equations` -- Is a "smaller" array with only the DOFs numbers of the prescribed equations.
    pub fn prescribed(&self, essential: &Essential) -> Result<(Vec<bool>, Vec<usize>), StrError> {
        let mut prescribed = vec![false; self.equations.n_equation];
        let mut p_equations = Vec::new();
        for (point_id, dof) in essential.all.keys() {
            let eq = self.equations.eq(*point_id, *dof)?;
            prescribed[eq] = true;
            p_equations.push(eq);
        }
        Ok((prescribed, p_equations))
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Data;
    use crate::base::{Ebc, Element, Essential, SampleParams};
    use gemlab::mesh::{Cell, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p2 = SampleParams::param_solid();
        assert_eq!(
            Data::new(&mesh, [(2, Element::Solid(p2))]).err(),
            Some("cannot find CellAttributeId in Attributes map")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        assert_eq!(data.equations.n_equation, 6);
    }

    #[test]
    fn n_local_eq_works() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        assert_eq!(data.n_local_eq(&mesh.cells[0]).unwrap(), 3);

        let wrong_cell = Cell {
            id: 0,
            attribute_id: 1,
            kind: GeoKind::Qua4,
            points: vec![0, 1, 2, 3],
        };
        assert_eq!(
            data.n_local_eq(&wrong_cell).err(),
            Some("cannot find (CellAttributeId, GeoKind) in ElementInfoMap")
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
        let data = Data::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let mut essential = Essential::new();
        let zero = |_| 0.0;
        assert_eq!(zero(1.0), 0.0);
        essential.at(&[0, 4], Ebc::Pl(zero));
        let (prescribed, mut p_equations) = data.prescribed(&essential).unwrap();
        assert_eq!(prescribed, &[true, false, false, false, true]);
        p_equations.sort();
        assert_eq!(p_equations, &[0, 4]);

        essential.at(&[3], Ebc::T(zero));
        assert_eq!(
            data.prescribed(&essential).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
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
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .at(&[0], Ebc::Ux(zero))
            .at(&[0], Ebc::Uy(zero))
            .at(&[1, 2], Ebc::Uy(zero));
        let (prescribed, mut p_equations) = data.prescribed(&essential).unwrap();
        assert_eq!(
            prescribed,
            //   0     1      2     3      4     5      6      7      8      9
            &[true, true, false, true, false, true, false, false, false, false]
        );
        p_equations.sort();
        assert_eq!(p_equations, &[0, 1, 3, 5]);
    }
}
