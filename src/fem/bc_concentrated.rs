use crate::base::{BcNatural, Schema};
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
    ///
    /// The function is `(stage, t) -> load`
    function: Option<&'a Box<dyn Fn(usize, f64) -> f64 + 'a>>,
}

/// Implements an array of BcConcentrated
pub struct BcConcentratedArray<'a> {
    /// All values
    pub all: Vec<BcConcentrated<'a>>,
}

impl<'a> BcConcentratedArray<'a> {
    /// Allocates a new instance
    pub fn new(schema: &Schema, natural: &'a BcNatural) -> Result<Self, StrError> {
        let mut all = Vec::with_capacity(natural.at_points.len() + 1);
        for (point_id, pbc, value, f_index) in &natural.at_points {
            let eq = schema.get_eq(*point_id, pbc.dof())?;
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

    /// Adds all concentrated load values at given time to the external forces vector F
    pub fn add_to_ff(&self, ff: &mut Vector, stage: usize, time: f64) {
        for e in &self.all {
            let value = match e.function {
                Some(f) => (f)(stage, time),
                None => e.value,
            };
            ff[e.eq] += value;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcConcentratedArray;
    use crate::base::{BcNatural, ParamSolid, Pbc, Schema};
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();

        let mut natural = BcNatural::new();
        natural.points(&[100], Pbc::Fx, -10.0);
        assert_eq!(
            BcConcentratedArray::new(&schema, &natural).err(),
            Some("cannot get equation number because point_id is out of bounds")
        );
    }

    #[test]
    fn add_to_ff_works() {
        let mesh = Samples::one_tet4();
        let p1 = ParamSolid::sample_linear_elastic();
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let mut natural = BcNatural::new();
        natural.points(&[0], Pbc::Fx, -20.0);
        natural.points(&[1], Pbc::Fy, -20.0);
        natural.points(&[2], Pbc::Fz, -20.0);
        let b_points = BcConcentratedArray::new(&schema, &natural).unwrap();
        let mut ff = Vector::new(4 * 3);
        b_points.add_to_ff(&mut ff, 0, 0.0);
        assert_eq!(
            ff.as_data(),
            &[
                -20.0, 0.0, 0.0, // 0
                0.0, -20.0, 0.0, // 1
                0.0, 0.0, -20.0, // 2
                0.0, 0.0, 0.0, // 3
            ]
        );
    }
}
