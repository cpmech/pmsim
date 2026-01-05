use super::{BcConcentratedArray, BcDistributedArray, BcPrescribed, Elements, FemState, LinearSystem};
use crate::base::{BcEssential, BcNatural, Config, Schema};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{vec_copy, vec_inner, vec_minus, Stopwatch};

/// Implements common (shared) functionality for all FEM solvers
pub(crate) struct SolverCommon<'a> {
    /// Holds the configuration
    pub(crate) config: &'a Config<'a>,

    /// Holds the material parameters, element attributes, and equation numbers
    pub(crate) base: &'a Schema,

    // Holds a collection of concentrated loads
    pub(crate) bc_concentrated: BcConcentratedArray<'a>,

    // Holds a collection of boundary integration data
    pub(crate) bc_distributed: BcDistributedArray<'a>,

    /// Holds a collection of prescribed (primary) values
    pub(crate) bc_prescribed: BcPrescribed<'a>,

    /// Holds a collection of elements
    pub(crate) elements: Elements<'a>,

    /// Holds variables to solve the global linear system
    pub(crate) ls: LinearSystem<'a>,

    /// Array to ignore prescribed equations when building the reduced system
    pub(crate) ignored_eqs: Vec<bool>,

    /// Unknown equation numbers
    pub(crate) unknown_eqs: Vec<usize>,

    /// Stopwatch to measure computer time
    pub(crate) stopwatch: Stopwatch,
}

impl<'a> SolverCommon<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        base: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
    ) -> Result<Self, StrError> {
        // check
        if let Some(msg) = config.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot allocate simulation because config.validate() failed");
        }

        // allocate auxiliary instances
        let bc_concentrated = BcConcentratedArray::new(base, natural)?;
        let bc_distributed = BcDistributedArray::new(mesh, base, config, natural)?;
        let bc_prescribed = BcPrescribed::new(base, essential)?;
        let elements = Elements::new(mesh, base, config)?;
        let linear_system = LinearSystem::new(base, config, &bc_prescribed, &elements, &bc_distributed)?;

        // array to ignore prescribed equations when building the reduced system
        let ndof = bc_prescribed.flags.len(); // number of DOFs = n_equation without Lagrange multipliers
        let ignored_eqs = if config.lagrange_mult_method {
            vec![false; ndof]
        } else {
            bc_prescribed.flags.clone()
        };

        // collect the unknown equations
        let neq_total = linear_system.neq_total;
        let unknown_eqs: Vec<_> = (0..neq_total)
            .filter(|&eq| config.lagrange_mult_method || !ignored_eqs[eq])
            .collect();

        // return new instance
        Ok(SolverCommon {
            config,
            base,
            bc_concentrated,
            bc_distributed,
            bc_prescribed,
            elements,
            ls: linear_system,
            ignored_eqs,
            unknown_eqs,
            stopwatch: Stopwatch::new(),
        })
    }

    /// Calculates Y (internal forces)
    pub fn calc_yy(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // clear vector
        self.ls.yy.fill(0.0);

        // calculate all element local vectors
        self.elements.assemble_yy(&mut self.ls.yy, state, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.bc_distributed
            .assemble_yy(&mut self.ls.yy, state, &self.ignored_eqs)?;
        Ok(())
    }

    /// Calculates F and ΔF
    ///
    /// Returns the load reversal flag
    ///
    /// ```text
    /// F_old := F(t)
    /// ΔF = F(t+Δt) - F(t)
    /// ```
    pub fn calc_ff_and_ddff(&mut self, step: usize, time: f64) -> Result<bool, StrError> {
        // make a copy of F and ΔF
        vec_copy(&mut self.ls.ff_old, &self.ls.ff).unwrap();
        vec_copy(&mut self.ls.ddff_old, &self.ls.ddff).unwrap();

        // update F ---------------------------------------------------------------

        // clear vector
        self.ls.ff.fill(0.0);

        // calculate all element local vectors
        self.elements
            .assemble_ff(&mut self.ls.ff, step, time, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.bc_distributed
            .assemble_ff(&mut self.ls.ff, step, time, &self.ignored_eqs)?;

        // add concentrated loads
        self.bc_concentrated.add_to_ff(&mut self.ls.ff, step, time);

        // ------------------------------------------------------------------------

        // calculate ΔF = F - F_old
        vec_minus(&mut self.ls.ddff, &self.ls.ff, &self.ls.ff_old).unwrap();

        // check if load reversal occurred
        let dot = vec_inner(&self.ls.ddff_old, &self.ls.ddff);
        let reverse = dot < 0.0 && self.config.consider_load_reversal;
        Ok(reverse)
    }

    /// Assembles the (augmented) global matrix K
    pub fn assemble_kk(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // reset pointer in K matrix == clear all values
        self.ls.kk.reset();

        // calculates all Ke matrices (local Jacobian matrix; derivative of Ye w.r.t u) and adds them to K
        self.elements.assemble_kk(&mut self.ls.kk, state, &self.ignored_eqs)?;
        self.bc_distributed
            .assemble_kk(&mut self.ls.kk, state, &self.ignored_eqs)?;
        Ok(())
    }

    /// Updates the (augmented) vectors of primary variables U, V, A
    pub fn update_primary_variables(&mut self, state: &mut FemState) -> Result<(), StrError> {
        let mdu = &mut self.ls.mdu;
        if self.config.transient {
            // update U, V, and ΔU vectors
            for i in &self.unknown_eqs {
                state.u[*i] -= mdu[*i];
                state.v[*i] = state.beta1 * state.u[*i] - state.u_star[*i];
                state.ddu[*i] -= mdu[*i];
            }
        } else {
            // update U and ΔU vectors
            for i in &self.unknown_eqs {
                state.u[*i] -= mdu[*i];
                state.ddu[*i] -= mdu[*i];
            }
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::SolverCommon;
    use crate::base::{BcEssential, BcNatural, Config, ParamSolid, Schema};
    use gemlab::mesh::Samples;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let mut schema = Schema::new();
        schema.add_solid(1, p1).build(&mesh).unwrap();
        let essential = BcEssential::new();
        let natural = BcNatural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt_min(-1.0);
        assert_eq!(
            SolverCommon::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
    }
}
