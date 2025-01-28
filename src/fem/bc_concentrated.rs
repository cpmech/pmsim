use super::FemInput;
use crate::base::{Natural, Pbc};
use crate::StrError;
use gemlab::mesh::PointId;
use russell_lab::Vector;

/// Assists in calculating concentrated loads
pub struct ConcentratedLoad {
    /// Point boundary condition
    pub pbc: Pbc,

    /// Equation corresponding to the concentrated load
    pub eq: usize,
}

/// Holds a collection of concentrated loads
pub struct ConcentratedLoads {
    /// All values
    pub all: Vec<ConcentratedLoad>,
}

impl ConcentratedLoad {
    /// Allocates new instance
    pub fn new(input: &FemInput, point_id: PointId, pbc: Pbc) -> Result<Self, StrError> {
        Ok(ConcentratedLoad {
            pbc,
            eq: input.equations.eq(point_id, pbc.dof())?,
        })
    }

    /// Adds the concentrated load value at given time to the global residual
    pub fn add_to_residual(&self, residual: &mut Vector, _time: f64) {
        let value = match self.pbc {
            Pbc::Fx(fx) => fx,
            Pbc::Fy(fy) => fy,
            Pbc::Fz(fz) => fz,
        };
        // note the negative sign
        residual[self.eq] -= value;
    }
}

impl ConcentratedLoads {
    /// Allocates new instance
    pub fn new(input: &FemInput, natural: &Natural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.at_points.len() + 1);
        for (point_id, pbc) in &natural.at_points {
            all.push(ConcentratedLoad::new(input, *point_id, *pbc)?);
        }
        Ok(ConcentratedLoads { all })
    }

    /// Adds all concentrated load values at given time to the global residual
    pub fn add_to_residual(&self, residual: &mut Vector, time: f64) {
        self.all.iter().for_each(|e| e.add_to_residual(residual, time));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ConcentratedLoad, ConcentratedLoads};
    use crate::base::{Etype, Natural, ParamSolid, Pbc};
    use crate::fem::FemInput;
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))]).unwrap();
        assert_eq!(
            ConcentratedLoad::new(&input, 123, Pbc::Fy(-10.0)).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );

        let mut natural = Natural::new();
        natural.points(&[100], Pbc::Fx(-10.0));
        assert_eq!(
            ConcentratedLoads::new(&input, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
    }

    #[test]
    fn add_to_residual_works() {
        let mesh = Samples::one_tet4();
        let p1 = ParamSolid::sample_linear_elastic();
        let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))]).unwrap();
        let mut natural = Natural::new();
        natural.points(&[0], Pbc::Fx(-20.0));
        natural.points(&[1], Pbc::Fy(-20.0));
        natural.points(&[2], Pbc::Fz(-20.0));
        let b_points = ConcentratedLoads::new(&input, &natural).unwrap();
        let mut residual = Vector::new(4 * 3);
        b_points.add_to_residual(&mut residual, 0.0);
        assert_eq!(
            residual.as_data(),
            &[
                20.0, 0.0, 0.0, // 0
                0.0, 20.0, 0.0, // 1
                0.0, 0.0, 20.0, // 2
                0.0, 0.0, 0.0, // 3
            ]
        );
    }
}
