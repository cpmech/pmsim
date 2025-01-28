use super::FemMesh;
use crate::base::Natural;
use crate::StrError;
use russell_lab::Vector;

/// Assists in calculating a concentrated load boundary condition
///
/// This data structure corresponds to a single Natural (Neumann) boundary condition
pub struct BcConcentrated<'a> {
    /// Equation corresponding to the concentrated load
    eq: usize,

    /// Specified BC value (overridden by the function, if not None)
    value: f64,

    /// Function to calculate the BC value (overrides the value, if not None)
    function: Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>,
}

/// Implements an array of BcConcentrated
pub struct BcConcentratedArray<'a> {
    /// All values
    pub all: Vec<BcConcentrated<'a>>,
}

impl<'a> BcConcentrated<'a> {
    /// Adds the concentrated load value at given time to the global residual
    pub fn add_to_residual(&self, residual: &mut Vector, time: f64) {
        let value = match self.function {
            Some(f) => (f)(time),
            None => self.value,
        };
        // note the negative sign
        residual[self.eq] -= value;
    }
}

impl<'a> BcConcentratedArray<'a> {
    /// Allocates a new instance
    pub fn new(fem: &FemMesh, natural: &'a Natural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.at_points.len() + 1);
        for (point_id, pbc, value, f_index) in &natural.at_points {
            let eq = fem.equations.eq(*point_id, pbc.dof())?;
            let function = match f_index {
                Some(index) => Some(&natural.functions[*index]),
                None => None,
            };
            all.push(BcConcentrated {
                eq,
                value: *value,
                function,
            });
        }
        Ok(BcConcentratedArray { all })
    }

    /// Adds all concentrated load values at given time to the global residual
    pub fn add_to_residual(&self, residual: &mut Vector, time: f64) {
        self.all.iter().for_each(|e| e.add_to_residual(residual, time));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcConcentratedArray;
    use crate::base::{Elem, Natural, ParamSolid, Pbc};
    use crate::fem::FemMesh;
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

        let mut natural = Natural::new();
        natural.points(&[100], Pbc::Fx, -10.0);
        assert_eq!(
            BcConcentratedArray::new(&fem, &natural).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
    }

    #[test]
    fn add_to_residual_works() {
        let mesh = Samples::one_tet4();
        let p1 = ParamSolid::sample_linear_elastic();
        let fem = FemMesh::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut natural = Natural::new();
        natural.points(&[0], Pbc::Fx, -20.0);
        natural.points(&[1], Pbc::Fy, -20.0);
        natural.points(&[2], Pbc::Fz, -20.0);
        let b_points = BcConcentratedArray::new(&fem, &natural).unwrap();
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
