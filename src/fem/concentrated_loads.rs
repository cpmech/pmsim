use super::Data;
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
    pub fn new(data: &Data, point_id: PointId, pbc: Pbc) -> Result<Self, StrError> {
        Ok(ConcentratedLoad {
            pbc,
            eq: data.equations.eq(point_id, pbc.dof())?,
        })
    }

    /// Adds the concentrated load value at given time to the global residual
    pub fn add_to_residual(&self, residual: &mut Vector, time: f64) {
        let value = match self.pbc {
            Pbc::Fx(f) => f(time),
            Pbc::Fy(f) => f(time),
            Pbc::Fz(f) => f(time),
        };
        // note the negative sign
        residual[self.eq] -= value;
    }
}

impl ConcentratedLoads {
    /// Allocates new instance
    pub fn new(data: &Data, natural: &Natural) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = natural
            .concentrated
            .iter()
            .map(|(point_id, pbc)| ConcentratedLoad::new(data, *point_id, *pbc))
            .collect();
        match res {
            Ok(all) => Ok(ConcentratedLoads { all }),
            Err(e) => Err(e),
        }
    }

    /// Adds all concentrated load values at given time to the global residual
    #[inline]
    pub fn add_to_residual(&self, residual: &mut Vector, time: f64) {
        self.all.iter().for_each(|e| e.add_to_residual(residual, time));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{ConcentratedLoad, ConcentratedLoads};
    use crate::base::{Element, Natural, Pbc, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let minus_ten = |_| -10.0;
        assert_eq!(minus_ten(0.0), -10.0);
        assert_eq!(
            ConcentratedLoad::new(&data, 123, Pbc::Fy(minus_ten)).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );

        let mut natural = Natural::new();
        natural.at(&[100], Pbc::Fx(minus_ten));
        assert_eq!(
            ConcentratedLoads::new(&data, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
    }

    #[test]
    fn add_to_residual_works() {
        let mesh = Samples::one_tet4();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut natural = Natural::new();
        let f = |_| -20.0;
        assert_eq!(f(0.0), -20.0);
        natural.at(&[0], Pbc::Fx(f));
        natural.at(&[1], Pbc::Fy(f));
        natural.at(&[2], Pbc::Fz(f));
        let b_points = ConcentratedLoads::new(&data, &natural).unwrap();
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
