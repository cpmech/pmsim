use super::FemState;
use crate::base::{BcEssential, Dof, Schema};
use crate::StrError;
use gemlab::mesh::PointId;
use russell_lab::Vector;
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};
use std::collections::HashMap;

/// Calculates the values associated with the prescribed essential boundary conditions
pub(crate) struct BcPrescribed<'a> {
    essential: &'a BcEssential<'a>,

    /// Manages equation numbers (prescribed versus unknown)
    pub(crate) handler: EquationHandler,

    /// Holds the pairs (PointId, Dof) for the prescribed equations (only)
    ///
    /// len = n_prescribed
    pairs: Vec<(PointId, Dof)>,
}

impl<'a> BcPrescribed<'a> {
    /// Allocates a new instance
    pub fn new(schema: &Schema, essential: &'a BcEssential) -> Result<Self, StrError> {
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

        // Return the instance
        Ok(BcPrescribed {
            essential,
            handler,
            pairs,
        })
    }

    /// Returns the value of the prescribed DOF at given time
    ///
    /// # Panics
    ///
    /// This function will panic if the equation number is out of bounds.
    pub fn value(&self, eq: usize, time: f64) -> f64 {
        let pair = self.pairs[eq];
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
        let neq = self.handler.neq();
        for ip in 0..self.handler.np() {
            let i = self.handler.prescribed()[ip];
            let j = neq + ip;
            let lag = state.u[j];
            let val = self.value(ip, state.time);
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
        let neq = self.handler.neq();
        let sym = kk.get_info().3;
        match sym {
            Sym::YesLower => {
                for ip in 0..self.handler.np() {
                    let i = self.handler.prescribed()[ip];
                    let j = neq + ip;
                    kk.put(j, i, 1.0).unwrap(); // A
                }
            }
            Sym::YesUpper => {
                for ip in 0..self.handler.np() {
                    let i = self.handler.prescribed()[ip];
                    let j = neq + ip;
                    kk.put(i, j, 1.0).unwrap(); // Aᵀ
                }
            }
            Sym::YesFull | Sym::No => {
                for ip in 0..self.handler.np() {
                    let i = self.handler.prescribed()[ip];
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
        for eq in self.handler.prescribed() {
            kk.put(*eq, *eq, 1.0).unwrap();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {}
