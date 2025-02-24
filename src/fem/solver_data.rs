use super::{BcConcentratedArray, BcDistributedArray, BcPrescribed, FemState};
use super::{Elements, FemBase, LinearSystem};
use crate::base::{Config, Essential, Natural};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::vec_add;

/// Holds data for FEM solvers
pub(crate) struct SolverData<'a> {
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
    pub(crate) ignored: Vec<bool>,

    /// Unknown equation numbers
    pub(crate) unknown: Vec<usize>,
}

impl<'a> SolverData<'a> {
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
        let ignore = if config.lagrange_mult_method {
            vec![false; ndof]
        } else {
            bc_prescribed.flags.clone()
        };

        // collect the unknown equations
        let neq_total = linear_system.neq_total;
        let unknown_equations: Vec<_> = (0..neq_total)
            .filter(|&eq| config.lagrange_mult_method || !ignore[eq])
            .collect();

        // return new instance
        Ok(SolverData {
            config,
            bc_concentrated,
            bc_distributed,
            bc_prescribed,
            elements,
            ls: linear_system,
            ignored: ignore,
            unknown: unknown_equations,
        })
    }

    /// Assembles the internal forces vector (F_int)
    pub fn assemble_ff_int(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // clear F_int vector
        self.ls.ff_int.fill(0.0);

        // calculate all element local vectors and add them to the global vectors
        self.elements
            .assemble_f_int(&mut self.ls.ff_int, state, &self.ignored)?;

        // calculate all boundary elements local vectors and add them to the global vectors
        self.bc_distributed
            .assemble_f_int(&mut self.ls.ff_int, state, &self.ignored)?;
        Ok(())
    }

    /// Assembles the external forces vector (F_ext)
    pub fn assemble_ff_ext(&mut self, t: f64) -> Result<(), StrError> {
        // clear F_ext vector
        self.ls.ff_ext.fill(0.0);

        // calculate all element local vectors and add them to the global vectors
        self.elements.assemble_f_ext(&mut self.ls.ff_ext, t, &self.ignored)?;

        // calculate all boundary elements local vectors and add them to the global vectors
        self.bc_distributed
            .assemble_f_ext(&mut self.ls.ff_ext, t, &self.ignored)?;

        // add concentrated loads to the external forces vector
        self.bc_concentrated.add_to_ff_ext(&mut self.ls.ff_ext, t);
        Ok(())
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
        self.elements.assemble_kke(&mut self.ls.kk, state, &self.ignored)?;
        self.bc_distributed
            .assemble_kke(&mut self.ls.kk, state, &self.ignored)?;
        Ok(())
    }

    /// Updates the (augmented) vectors of primary variables U, V, A
    pub fn update_primary_variables(&mut self, state: &mut FemState) -> Result<(), StrError> {
        let mdu = &mut self.ls.mdu;
        if self.config.transient {
            // update U, V, and ΔU vectors
            for i in &self.unknown {
                state.u[*i] -= mdu[*i];
                state.v[*i] = state.beta1 * state.u[*i] - state.u_star[*i];
                state.ddu[*i] -= mdu[*i];
            }
        } else {
            // update U and ΔU vectors
            for i in &self.unknown {
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
    use super::SolverData;
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
        config.set_dt_min(-1.0);
        assert_eq!(
            SolverData::new(&mesh, &base, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
    }
}
