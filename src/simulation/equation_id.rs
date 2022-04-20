use super::{Configuration, Dof, NDOF_PER_NODE_TOTAL};
use crate::StrError;
use gemlab::mesh::PointId;
use russell_lab::NumMatrix;

/// Holds equation identification numbers (all DOF numbers)
pub struct EquationId {
    /// Total number of equations
    ///
    /// The maximum number of equations is 2³¹ - 1 = 2,147,483,647
    /// (~2.1 billion! which will break any linear solver)
    nequation: i32,

    /// Matrix of equation identification numbers (npoint,ndof)
    ///
    /// * numbers are stored here following the **one-based** convention; i.e., they starts from 1
    /// * `eid < 0`: negative numbers indicate that the point/DOF has a prescribed (imposed) value
    /// * `eid = 0`: zero means that the point/DOF is not being used (unassigned)
    /// * `eid > 0`: positive numbers indicate that the point/DOF corresponds to an unknown value
    eid_one_based: NumMatrix<i64>,

    /// Holds a subset of equation identification numbers corresponding to prescribed DOFs
    ///
    /// Note: this vector holds **zero-based** indices.
    prescribed: Vec<usize>,
}

impl EquationId {
    /*
    /// Allocates a new instance
    ///
    /// This function will also initialize the equation_id matrix with the
    /// equation identification numbers of prescribed (imposed) DOFs
    pub fn new(config: &Configuration) -> Self {
        // initialize the eid matrix with UNASSIGNED values
        let npoint = config.mesh.points.len();
        let mut eid_one_based = NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, 0);

        // sort the prescribed (point_id,dof) pairs to have a deterministic numbering
        let mut pairs: Vec<_> = config.essential_bcs.keys().collect();
        pairs.sort_by(|a, b| a.cmp(&b));

        // set the (one-based and negative) eid of prescribed equations
        let mut nequation: i32 = 0;
        let mut prescribed = vec![0; config.essential_bcs.len()];
        let mut i = 0;
        for (point_id, dof) in pairs {
            nequation += 1;
            eid_one_based[*point_id][*dof as usize] = -nequation as i64;
            prescribed[i] = nequation as usize - 1;
            i += 1;
        }

        // done
        EquationId {
            nequation,
            eid_one_based,
            prescribed,
        }
    }
    */
    pub fn new(npoint: usize) -> Self {
        EquationId {
            nequation: 0,
            eid_one_based: NumMatrix::filled(npoint, NDOF_PER_NODE_TOTAL, 0),
            prescribed: Vec::new(),
        }
    }

    /// Activates an equation corresponding to a point/DOF pair
    ///
    /// This function will also increment the number of equations,
    /// if the equation does not exist yet.
    ///
    /// # Output
    ///
    /// * `(eid,prescribed)` -- A pair with the equation identification number and
    ///                         whether the `eid` corresponds to a prescribed DOF or not.
    pub fn activate(&mut self, point_id: PointId, dof: Dof) -> (usize, bool) {
        let eid_one_based = self.eid_one_based[point_id][dof as usize];
        if eid_one_based < 0 {
            (-eid_one_based as usize - 1, true) // prescribed
        } else if eid_one_based == 0 {
            self.nequation += 1;
            self.eid_one_based[point_id][dof as usize] = self.nequation as i64;
            (self.nequation as usize - 1, false) // not prescribed; newly activated
        } else {
            (eid_one_based as usize - 1, false) // not prescribed; previously activated
        }
    }

    /// Returns the equation identification number corresponding to a point/DOF pair
    ///
    /// # Output
    ///
    /// * `(eid,prescribed)` -- A pair with the equation identification number and
    ///                         whether the `eid` corresponds to a prescribed DOF or not.
    pub fn eid(&self, point_id: PointId, dof: Dof) -> Result<(usize, bool), StrError> {
        let eid_one_based = self.eid_one_based[point_id][dof as usize];
        if eid_one_based < 0 {
            Ok((-eid_one_based as usize - 1, true)) // prescribed
        } else if eid_one_based == 0 {
            Err("(point_id,dof) pair does not have an assigned equation id yet")
        } else {
            Ok((eid_one_based as usize - 1, false)) // not prescribed; i.e., unknown DOF
        }
    }

    /// Accesses the equation identification numbers corresponding to prescribed DOFs
    pub fn prescribed(&self) -> &Vec<usize> {
        &self.prescribed
    }

    /// Returns the number of points
    pub fn npoint(&self) -> usize {
        self.eid_one_based.nrow()
    }

    /// Returns the current total number of equations (prescribed and unknown DOFs)
    pub fn nequation(&self) -> usize {
        self.nequation as usize
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::EquationId;
    use crate::simulation::{Configuration, Dof, NDOF_PER_NODE_TOTAL};
    use crate::StrError;
    use gemlab::mesh::Mesh;

    fn mesh_square() -> Mesh {
        Mesh::from_text(
            r"
            #  3------2
            #  |      |
            #  |      |
            #  0------1

            # space_ndim npoint ncell
                       2      4     1

            # id   x   y
               0 0.0 0.0
               1 1.0 0.0
               2 1.0 1.0
               3 0.0 1.0

            # id att geo_ndim nnode  point_ids...
               0   1        2     4  0 1 2 3",
        )
        .unwrap()
    }

    /*

    #[test]
    fn new_works() -> Result<(), StrError> {
        let mesh = mesh_square();
        let mut config = Configuration::new(&mesh);

        let eqs = EquationId::new(&config);
        // (0,0) 3------2 (0,0) (ux,uy)
        //       |      |       one-based equation
        //       |      |       identification numbers
        // (0,0) 0------1 (0,0) (negative means prescribed)
        assert_eq!(eqs.nequation, 0);
        assert_eq!(eqs.eid_one_based.dims(), (4, NDOF_PER_NODE_TOTAL));
        assert_eq!(eqs.prescribed.len(), 0);
        assert_eq!(
            format!("{}", eqs.eid_one_based),
            "┌                     ┐\n\
             │ 0 0 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 0 0 │\n\
             │ 0 0 0 0 0 0 0 0 0 0 │\n\
             └                     ┘"
        );

        config.ebc_points(&[0, 3], &[Dof::Ux, Dof::Uy], |_, _| 0.0)?;

        let eqs = EquationId::new(&config);
        // (-3,-4) 3------2 (0,0) (ux,uy)
        //         |      |       one-based equation
        //         |      |       identification numbers
        // (-1,-2) 0------1 (0,0) (negative means prescribed)
        assert_eq!(eqs.nequation, 4);
        assert_eq!(eqs.eid_one_based.dims(), (4, NDOF_PER_NODE_TOTAL));
        assert_eq!(eqs.prescribed, [0, 1, 2, 3]);
        assert_eq!(
            format!("{}", eqs.eid_one_based),
            "┌                               ┐\n\
             │ -1 -2  0  0  0  0  0  0  0  0 │\n\
             │  0  0  0  0  0  0  0  0  0  0 │\n\
             │  0  0  0  0  0  0  0  0  0  0 │\n\
             │ -3 -4  0  0  0  0  0  0  0  0 │\n\
             └                               ┘"
        );
        Ok(())
    }

    #[test]
    fn activate_works() -> Result<(), StrError> {
        let mesh = mesh_square();
        let mut config = Configuration::new(&mesh);
        config.ebc_points(&[0, 3], &[Dof::Ux, Dof::Uy], |_, _| 0.0)?;

        let mut eqs = EquationId::new(&config);
        // (-3,-4) 3------2 (0,0) (ux,uy)
        //         |      |       one-based equation
        //         |      |       identification numbers
        // (-1,-2) 0------1 (0,0) (negative means prescribed)
        assert_eq!(eqs.nequation, 4);
        assert_eq!(eqs.prescribed, [0, 1, 2, 3]);

        // existent (prescribed)
        let (eid, prescribed) = eqs.activate(0, Dof::Ux);
        assert_eq!(eid, 0); // 0 = -(-1) -1
        assert_eq!(prescribed, true);
        assert_eq!(eqs.nequation, 4);

        // new (not prescribed)
        let (eid, prescribed) = eqs.activate(1, Dof::Ux);
        // (-3,-4) 3------2 (0,0) (ux,uy)
        //         |      |       one-based equation
        //         |      |       identification numbers
        // (-1,-2) 0------1 (5,0) (negative means prescribed)
        assert_eq!(eid, 4); // 4 = 5 - 1
        assert_eq!(prescribed, false);
        assert_eq!(eqs.nequation, 5);

        // new (not prescribed)
        let (eid, prescribed) = eqs.activate(2, Dof::Uy);
        // (-3,-4) 3------2 (0,6) (ux,uy)
        //         |      |       one-based equation
        //         |      |       identification numbers
        // (-1,-2) 0------1 (5,0) (negative means prescribed)
        assert_eq!(eid, 5); // 5 = 6 - 1
        assert_eq!(prescribed, false);
        assert_eq!(eqs.nequation, 6);

        // existent (not prescribed)
        let (eid, prescribed) = eqs.activate(1, Dof::Ux);
        assert_eq!(eid, 4); // 4 = 5 - 1
        assert_eq!(prescribed, false);

        // check matrix
        assert_eq!(
            format!("{}", eqs.eid_one_based),
            "┌                               ┐\n\
             │ -1 -2  0  0  0  0  0  0  0  0 │\n\
             │  5  0  0  0  0  0  0  0  0  0 │\n\
             │  0  6  0  0  0  0  0  0  0  0 │\n\
             │ -3 -4  0  0  0  0  0  0  0  0 │\n\
             └                               ┘"
        );
        Ok(())
    }

    #[test]
    fn getters_work() -> Result<(), StrError> {
        let mesh = mesh_square();
        let mut config = Configuration::new(&mesh);
        config.ebc_points(&[0, 3], &[Dof::Ux, Dof::Uy], |_, _| 0.0)?;

        let mut eqs = EquationId::new(&config);
        eqs.activate(1, Dof::Ux);
        eqs.activate(2, Dof::Uy);
        // (-3,-4) 3------2 (0,6) (ux,uy)
        //         |      |       one-based equation
        //         |      |       identification numbers
        // (-1,-2) 0------1 (5,0) (negative means prescribed)
        assert_eq!(eqs.nequation, 6);
        assert_eq!(eqs.prescribed, [0, 1, 2, 3]);
        assert_eq!(
            format!("{}", eqs.eid_one_based),
            "┌                               ┐\n\
             │ -1 -2  0  0  0  0  0  0  0  0 │\n\
             │  5  0  0  0  0  0  0  0  0  0 │\n\
             │  0  6  0  0  0  0  0  0  0  0 │\n\
             │ -3 -4  0  0  0  0  0  0  0  0 │\n\
             └                               ┘"
        );

        // check
        assert_eq!(eqs.eid(0, Dof::Ux), Ok((0, true)));
        assert_eq!(eqs.eid(0, Dof::Uy), Ok((1, true)));
        assert_eq!(eqs.eid(3, Dof::Ux), Ok((2, true)));
        assert_eq!(eqs.eid(3, Dof::Uy), Ok((3, true)));
        assert_eq!(eqs.eid(1, Dof::Ux), Ok((4, false)));
        assert_eq!(eqs.eid(2, Dof::Uy), Ok((5, false)));
        assert_eq!(
            eqs.eid(1, Dof::Uy).err(),
            Some("(point_id,dof) pair does not have an assigned equation id yet")
        );
        assert_eq!(
            eqs.eid(2, Dof::Ux).err(),
            Some("(point_id,dof) pair does not have an assigned equation id yet")
        );
        assert_eq!(eqs.prescribed(), &[0, 1, 2, 3]);
        assert_eq!(eqs.npoint(), 4);
        assert_eq!(eqs.nequation(), 6);
        Ok(())
    }
    */
}
