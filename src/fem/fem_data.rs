use super::{BcConcentratedArray, BcDistributedArray, Elements, FemState, LinearSystem};
use crate::base::{BcEssential, BcNatural, Config, Dof, Schema};
use crate::StrError;
use gemlab::mesh::{Mesh, PointId};
use russell_lab::{vec_copy, vec_inner, vec_minus, Stopwatch, Vector};
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};
use std::collections::HashMap;

/// Implements common (shared) functionality for all FEM solvers
pub(crate) struct FemData<'a> {
    essential: &'a BcEssential<'a>,

    /// Holds the pairs (PointId, Dof) for the prescribed equations (only)
    ///
    /// len = n_prescribed
    presc_pairs: Vec<(PointId, Dof)>,

    /// Manages equation numbers (prescribed versus unknown)
    pub(crate) eq_handler: EquationHandler,

    /// Holds the configuration
    pub(crate) config: &'a Config<'a>,

    /// Holds element types, material parameters, and specifies the DOF numbering schema
    pub(crate) schema: &'a Schema,

    // Holds a collection of concentrated loads
    pub(crate) bc_concentrated: BcConcentratedArray<'a>,

    // Holds a collection of boundary integration data
    pub(crate) bc_distributed: BcDistributedArray<'a>,

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

impl<'a> FemData<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        essential: &'a BcEssential,
        natural: &'a BcNatural,
    ) -> Result<Self, StrError> {
        // check
        if let Some(msg) = config.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot allocate simulation because config.validate() failed");
        }

        // Collect the list of prescribed equations
        let n_prescribed = essential.size();
        let mut eq_to_pair = HashMap::new();
        let mut p_list = Vec::with_capacity(n_prescribed);
        for (point_id, dof) in essential.keys() {
            let eq = schema.get_eq(*point_id, *dof)?;
            p_list.push(eq);
            eq_to_pair.insert(eq, (*point_id, *dof));
        }

        // Allocate the equations handler
        let mut handler = EquationHandler::new(schema.get_neq()?);
        handler.recompute(&p_list);

        // Allocate array of (point_id, dof) pairs corresponding to prescribed equation numbers
        let mut pairs = Vec::with_capacity(n_prescribed);
        for eq in handler.prescribed() {
            pairs.push(eq_to_pair.get(eq).unwrap().to_owned());
        }

        // allocate auxiliary instances
        let bc_concentrated = BcConcentratedArray::new(schema, natural)?;
        let bc_distributed = BcDistributedArray::new(mesh, schema, config, natural)?;
        let elements = Elements::new(mesh, schema, config)?;
        let linear_system = LinearSystem::new(n_prescribed, schema, config, &elements, &bc_distributed)?;

        // array to ignore prescribed equations when building the reduced system
        let ndof = handler.neq(); // number of DOFs = n_equation without Lagrange multipliers
        let mut ignored_eqs = vec![false; ndof];
        if !config.lagrange_mult_method {
            for eq in handler.prescribed() {
                ignored_eqs[*eq] = true;
            }
        };

        // collect the unknown equations
        let neq_total = linear_system.neq_total;
        let unknown_eqs: Vec<_> = (0..neq_total)
            .filter(|&eq| config.lagrange_mult_method || !ignored_eqs[eq])
            .collect();

        // return new instance
        Ok(FemData {
            essential,
            presc_pairs: pairs,
            eq_handler: handler,
            config,
            schema,
            bc_concentrated,
            bc_distributed,
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

    /// Returns the value of the prescribed DOF at given time
    ///
    /// # Panics
    ///
    /// This function will panic if the equation number is out of bounds.
    pub fn prescribed_value(&self, eq: usize, time: f64) -> f64 {
        let pair = self.presc_pairs[eq];
        self.essential.value(pair.0, pair.1, time)
    }

    /// Assembles the contribution due to the prescribed DOFs into the global R vector (LMM)
    ///
    /// **LMM** means Lagrange Multiplier Method
    ///
    /// This function adds `Aᵀλ` to the global R vector at the non-prescribed equations and
    /// **sets** the prescribed equations to `A u - c`. Here, `c` is the prescribed value.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌         ┐
    ///  │  K   Aᵀ │ │ -δu │   │ R + Aᵀλ │
    ///  │         │ │     │ = │         │
    ///  │  A   0  │ │ -δλ │   │ A u - c │
    ///  └         ┘ └     ┘   └         ┘
    /// ```
    pub fn assemble_rr_lmm(&self, rr: &mut Vector, state: &FemState) {
        let neq = self.eq_handler.neq();
        for ip in 0..self.eq_handler.np() {
            let i = self.eq_handler.prescribed()[ip];
            let j = neq + ip;
            let lag = state.u[j];
            let val = self.prescribed_value(ip, state.time);
            rr[i] += lag; // Aᵀ λ  →  1 * λ
            rr[j] = state.u[i] - val; // A u - c  →  1 * u - c
        }
    }

    /// Assembles the constraint matrix into the global K matrix (LMM)
    ///
    /// **LMM** means Lagrange Multiplier Method
    ///
    /// This function adds the constraints matrix (Aᵀ and A) to K.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌         ┐
    ///  │  K   Aᵀ │ │ -δu │   │ R + Aᵀλ │
    ///  │         │ │     │ = │         │
    ///  │  A   0  │ │ -δλ │   │ A u - c │
    ///  └         ┘ └     ┘   └         ┘
    /// ```
    pub fn assemble_kk_lmm(&self, kk: &mut CooMatrix) {
        let neq = self.eq_handler.neq();
        let sym = kk.get_info().3;
        match sym {
            Sym::YesLower => {
                for ip in 0..self.eq_handler.np() {
                    let i = self.eq_handler.prescribed()[ip];
                    let j = neq + ip;
                    kk.put(j, i, 1.0).unwrap(); // A
                }
            }
            Sym::YesUpper => {
                for ip in 0..self.eq_handler.np() {
                    let i = self.eq_handler.prescribed()[ip];
                    let j = neq + ip;
                    kk.put(i, j, 1.0).unwrap(); // Aᵀ
                }
            }
            Sym::YesFull | Sym::No => {
                for ip in 0..self.eq_handler.np() {
                    let i = self.eq_handler.prescribed()[ip];
                    let j = neq + ip;
                    kk.put(i, j, 1.0).unwrap(); // Aᵀ
                    kk.put(j, i, 1.0).unwrap(); // A
                }
            }
        }
    }

    /// Updates the diagonal of the global K matrix (RSM)
    ///
    /// **RSM** means Reduced-System Method
    ///
    /// This function put ones on the diagonal entries corresponding to the prescribed DOFs.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌   ┐
    ///  │  K   0  │ │ -δu │   │ R │
    ///  │         │ │     │ = │   │
    ///  │  0   1  │ │  0  │   │ 0 │
    ///  └         ┘ └     ┘   └   ┘
    /// ```
    ///
    /// Note that the prescribed values are zero (homogeneous BCs).
    pub fn assemble_kk_rsm(&self, kk: &mut CooMatrix) {
        for eq in self.eq_handler.prescribed() {
            kk.put(*eq, *eq, 1.0).unwrap();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemData;
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
            FemData::new(&mesh, &schema, &config, &essential, &natural).err(),
            Some("cannot allocate simulation because config.validate() failed")
        );
    }
}
