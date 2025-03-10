use super::{BcConcentratedArray, BcDistributedArray, BcPrescribed, FemState};
use super::{Elements, FemBase, LinearSystem};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::{vec_add, vec_copy, vec_inner, vec_minus, vec_update, Stopwatch};

/// Implements common (shared) functionality for all FEM solvers
pub(crate) struct SolverCommon<'a> {
    /// Holds the configuration
    config: &'a Config<'a>,

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
        base: &'a FemBase,
        config: &'a Config,
        essential: &'a Essential,
        natural: &'a Natural,
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

        // show information
        if config.verbose_timesteps || config.verbose_iterations {
            println!("\nINFORMATION ================================================================");
            println!("\n{}", linear_system.get_info());
        }

        // return new instance
        Ok(SolverCommon {
            config,
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

    /// Assembles the internal forces vector (F_int)
    pub fn assemble_ff_int(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // clear F_int vector
        self.ls.ff_int.fill(0.0);

        // calculate all element local vectors and add them to F_int
        self.elements
            .assemble_f_int(&mut self.ls.ff_int, state, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors and add them to F_int
        self.bc_distributed
            .assemble_f_int(&mut self.ls.ff_int, state, &self.ignored_eqs)?;
        Ok(())
    }

    /// Assembles the external forces vector (F_ext)
    ///
    /// Returns the load reversal flag
    ///
    /// ```text
    ///         ⎧ F_ext_old + λ ΔF_ext  if quasi-static/steady
    /// F_ext = ⎨
    ///         ⎩ F_ext(t)              if transient/dynamics
    /// ```
    pub fn assemble_ff_ext(&mut self, stage: usize, lambda: f64, t: f64) -> Result<bool, StrError> {
        let reverse = if self.config.steady {
            // make a copy of ΔF_ext
            vec_copy(&mut self.ls.ddff_ext_old, &self.ls.ddff_ext).unwrap();

            // assemble F_ext into tmp ------------------------------------------------

            // clear tmp vector
            self.ls.tmp.fill(0.0);

            // calculate all element local vectors and add them to tmp
            self.elements.assemble_f_ext(&mut self.ls.tmp, t, &self.ignored_eqs)?;

            // calculate all boundary elements local vectors and add them to tmp
            self.bc_distributed
                .assemble_f_ext(&mut self.ls.tmp, stage, t, &self.ignored_eqs)?;

            // add concentrated loads to tmp
            self.bc_concentrated.add_to_ff_ext(&mut self.ls.tmp, stage, t);

            // ------------------------------------------------------------------------

            // calculate ΔF_ext = tmp - F_ext
            vec_minus(&mut self.ls.ddff_ext, &self.ls.tmp, &self.ls.ff_ext).unwrap();

            // calculate F_ext += λ ΔF_ext
            vec_update(&mut self.ls.ff_ext, lambda, &self.ls.ddff_ext).unwrap();

            // check if load reversal occurred
            let dot = vec_inner(&self.ls.ddff_ext_old, &self.ls.ddff_ext);
            dot < 0.0 && self.config.consider_load_reversal
        } else {
            // assemble F_ext ---------------------------------------------------------

            // clear F_ext vector
            self.ls.ff_ext.fill(0.0);

            // calculate all element local vectors and add them to F_ext
            self.elements
                .assemble_f_ext(&mut self.ls.ff_ext, t, &self.ignored_eqs)?;

            // calculate all boundary elements local vectors and add them to F_ext
            self.bc_distributed
                .assemble_f_ext(&mut self.ls.ff_ext, stage, t, &self.ignored_eqs)?;

            // add concentrated loads to F_ext
            self.bc_concentrated.add_to_ff_ext(&mut self.ls.ff_ext, stage, t);

            // ------------------------------------------------------------------------

            // ignore load reversal
            false
        };
        Ok(reverse)
    }

    /// Calculates the residual vector R
    ///
    /// ```text
    /// R = F_int - lf * F_ext
    /// ```
    ///
    /// where `lf` is the loading factor.
    pub fn calculate_residuals_vector(&mut self, loading_factor: f64) {
        // R = F_int - lf * F_ext
        vec_add(&mut self.ls.rr, 1.0, &self.ls.ff_int, -loading_factor, &self.ls.ff_ext).unwrap();
    }

    /// Assembles the (augmented) global matrix K
    pub fn assemble_kk(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // reset pointer in K matrix == clear all values
        self.ls.kk.reset();

        // calculates all Ke matrices (local Jacobian matrix; derivative of f_int w.r.t u) and adds them to K
        self.elements.assemble_kke(&mut self.ls.kk, state, &self.ignored_eqs)?;
        self.bc_distributed
            .assemble_kke(&mut self.ls.kk, state, &self.ignored_eqs)?;
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
    use crate::base::{Config, Elem, Essential, Natural, ParamSolid};
    use crate::fem::FemBase;
    use gemlab::mesh::Samples;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_hex8();
        let mut p1 = ParamSolid::sample_linear_elastic();
        p1.ngauss = Some(123); // wrong
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let essential = Essential::new();
        let natural = Natural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_transient().set_ddt_min(-1.0);
        assert_eq!(
            SolverCommon::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
    }
}
