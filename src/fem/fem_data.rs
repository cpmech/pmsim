#![allow(unused)]

use super::{ElementsBoundary, ElementsInterior, FemState, LinearSystem, OutputFiles};
use crate::base::{BcEssential, BcNatural, Config, Dof, Schema};
use crate::StrError;
use gemlab::mesh::{Mesh, PointId};
use russell_lab::{vec_copy, vec_inner, vec_minus, Stopwatch, Vector};
use russell_nonlin::Stop;
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};
use std::collections::HashMap;
use std::fmt::Write;
use std::sync::Arc;
use uuid::Uuid;

/// Implements common (shared) functionality for all FEM solvers
pub struct FemData<'a> {
    /// Holds element types, material parameters, and specifies the DOF numbering schema
    pub(crate) schema: &'a Schema,

    /// Holds the configuration
    pub(crate) config: &'a Config<'a>,

    /// Manages equation numbers (prescribed versus unknown)
    pub(crate) eq_handler: EquationHandler,

    /// Holds the functions to calculate the prescribed values
    ///
    /// len = n_prescribed; use eq_handler.ip() to access an entry in this array
    pub(crate) presc_values: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>>,

    /// Holds pairs of (eq, fn) to calculate concentrated loads
    pub(crate) conc_loads: Vec<(usize, Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>)>,

    // Holds a collection of boundary elements
    pub(crate) boundaries: ElementsBoundary<'a>,

    /// Holds a collection of elements
    pub(crate) elements: ElementsInterior<'a>,

    /// Holds variables to solve the global linear system
    pub(crate) ls: LinearSystem<'a>,

    /// Array to ignore prescribed equations when building the reduced system
    pub(crate) ignored_eqs: Vec<bool>,

    /// Unknown equation numbers
    pub(crate) unknown_eqs: Vec<usize>,

    /// Handles output files
    pub(crate) files: OutputFiles,

    /// Stopwatch to measure computer time
    pub(crate) stopwatch: Stopwatch,

    /// Vector of internal forces
    ///
    /// dim = neq = nu + np
    pub(crate) yy: Vector,

    /// Vector of external forces
    ///
    /// dim = neq = nu + np
    pub(crate) ff: Vector,

    pub(crate) uuid: Uuid,
    pub(crate) state: FemState,
    pub(crate) neq: usize,
    pub(crate) nu: usize,
    pub(crate) np: usize,
    pub(crate) ndim: usize,
    pub(crate) sym: Sym,
    pub(crate) nnz_kk: usize,
    pub(crate) nnz_kk_bar: usize,
    pub(crate) nnz_kk_check: usize,
    pub(crate) kk_check: CooMatrix,
    pub(crate) u_check: Vector,
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
        // Check
        if let Some(msg) = config.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot start simulation because config.validate() failed");
        }

        // Start stopwatch
        let mut stopwatch = Stopwatch::new();

        // Generate the list of prescribed equations and a map from equation to (PointId, Dof)
        let n_prescribed = essential.functions.len();
        let mut p_list = Vec::with_capacity(n_prescribed);
        let mut eq_to_dof = HashMap::with_capacity(n_prescribed);
        for (point_id, dof) in essential.functions.keys() {
            let eq = schema.get_eq(*point_id, *dof)?;
            p_list.push(eq);
            eq_to_dof.insert(eq, (*point_id, *dof));
        }

        // Allocate the equations handler
        let mut eq_handler = EquationHandler::new(schema.get_neq()?);
        eq_handler.recompute(&p_list);

        // Allocate array of functions to calculate prescribed values
        let mut presc_values = Vec::with_capacity(n_prescribed);
        for eq in eq_handler.prescribed() {
            let point_dof = eq_to_dof.get(eq).unwrap();
            let f = essential.functions.get(point_dof).unwrap();
            presc_values.push(f.clone());
        }

        // Allocate array of concentrated loads
        let mut conc_loads = Vec::with_capacity(natural.at_points.len());
        for (point_id, pbc, f) in &natural.at_points {
            let eq = schema.get_eq(*point_id, pbc.dof())?;
            conc_loads.push((eq, f.clone()));
        }

        // Allocate auxiliary instances
        let boundaries = ElementsBoundary::new(mesh, schema, config, natural)?;
        let mut elements = ElementsInterior::new(mesh, schema, config)?;
        let linear_system = LinearSystem::new(n_prescribed, schema, config, &elements, &boundaries)?;

        // Array to ignore prescribed equations when building the reduced system
        let neq = eq_handler.neq(); // number of DOFs (without Lagrange multipliers)
        let mut ignored_eqs = vec![false; neq];
        if !config.lagrange_mult_method {
            for eq in eq_handler.prescribed() {
                ignored_eqs[*eq] = true;
            }
        };

        // Collect the unknown equations
        let neq_total = linear_system.neq_total;
        let unknown_eqs: Vec<_> = (0..neq_total)
            .filter(|&eq| config.lagrange_mult_method || !ignored_eqs[eq])
            .collect();

        // Allocate output files handler
        let mut files = OutputFiles::new(mesh, schema, config)?;

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        let mut state = FemState::new(&mesh, &schema, &essential, &config)?;

        // Initialize internal variables
        elements.initialize_internal_values(&mut state)?;

        // First output (must occur after initialize_internal_values)
        files.write_state(&config, &state)?;
        files.save_selected(&config, &schema, &state)?;

        // Determine if the global stiffness matrix is symmetric and it's enabled
        let symmetric = !config.ignore_symmetry && elements.all_sym_kk() && boundaries.all_sym_kk();

        // Determine symmetry type of the global stiffness matrix
        let genie = config.lin_sol_genie;
        let sym = genie.get_sym(symmetric);

        // Determine the system dimension
        let neq = eq_handler.neq();
        let nu = eq_handler.nu();
        let np = eq_handler.np();
        let ndim = if config.lagrange_mult_method { neq + np } else { nu };

        // Calculate the number of non-zero entries in the global stiffness matrix
        let mut nnz_kk = 0;
        let mut nnz_kk_bar = 0;
        let mut nnz_kk_check = 0;
        if config.lagrange_mult_method {
            elements.add_nnz_lmm(&mut nnz_kk, sym);
            boundaries.add_nnz_lmm(&mut nnz_kk, sym);
            if sym.triangular() {
                nnz_kk += np;
            } else {
                nnz_kk += 2 * np;
            }
        } else {
            elements.add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &eq_handler);
            boundaries.add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &eq_handler);
        }

        // Allocate K-check matrix for SPS
        let kk_check = if config.lagrange_mult_method || n_prescribed == 0 {
            CooMatrix::new(1, 1, 1, Sym::No).unwrap() // empty
        } else {
            CooMatrix::new(nu, np, nnz_kk_check, Sym::No).unwrap()
        };

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        // return new instance
        Ok(FemData {
            schema,
            config,
            eq_handler,
            presc_values,
            conc_loads,
            boundaries,
            elements,
            ls: linear_system,
            ignored_eqs,
            unknown_eqs,
            files,
            stopwatch,
            yy: Vector::new(neq),
            ff: Vector::new(neq),
            //
            uuid: Uuid::new_v4(),
            state,
            neq,
            nu,
            np,
            ndim,
            sym,
            nnz_kk,
            nnz_kk_bar,
            nnz_kk_check,
            kk_check,
            u_check: Vector::new(np),
        })
    }

    pub fn get_u_index(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        let eq = self.schema.get_eq(point_id, dof)?;
        if self.config.lagrange_mult_method {
            Ok(eq)
        } else {
            Ok(self.eq_handler.iu(eq))
        }
    }

    pub fn get_state(&self) -> &FemState {
        &self.state
    }

    pub fn write_state(&mut self) -> Result<(), StrError> {
        self.files.write_state(self.config, &self.state)
    }

    /// Prints information about the system
    pub fn print_system_info(&self, continuation: &str) {
        if self.config.verbose {
            let mut b = vec![vec![String::new(); 3]; 3];
            write!(&mut b[0][0], "neq  = {:?}", self.neq).unwrap();
            write!(&mut b[1][0], "np   = {:?}", self.np).unwrap();
            write!(&mut b[2][0], "ndim = {:?}", self.ndim).unwrap();
            write!(&mut b[0][1], "nnz(K)     = {:?}", self.nnz_kk).unwrap();
            write!(&mut b[1][1], "nnz(K-bar) = {:?}", self.nnz_kk_bar).unwrap();
            write!(&mut b[2][1], "sym(K)     = {:?}", self.sym).unwrap();
            write!(&mut b[0][2], "genie        = {:?}", self.config.lin_sol_genie).unwrap();
            write!(&mut b[1][2], "continuation = {}", continuation).unwrap();
            write!(
                &mut b[2][2],
                "EBC handler  = {}",
                if self.config.lagrange_mult_method { "LMM" } else { "SPS" }
            )
            .unwrap();
            let mut w = vec![0; 3];
            for i in 0..3 {
                for j in 0..3 {
                    w[j] = usize::max(w[j], b[i][j].len());
                }
            }
            let mut buf = String::new();
            for i in 0..3 {
                if i > 0 {
                    write!(&mut buf, "\n").unwrap();
                }
                for j in 0..3 {
                    if j > 0 {
                        write!(&mut buf, " │ ").unwrap();
                    }
                    write!(&mut buf, "{:1$}", b[i][j], w[j]).unwrap();
                }
            }
            write!(&mut buf, "\n").unwrap();
            println!("\n{}", buf);
        }
    }

    /// Calculates Y (internal forces)
    pub fn calc_yy(&mut self) -> Result<(), StrError> {
        // clear vector
        self.yy.fill(0.0);

        // calculate all element local vectors
        self.elements
            .assemble_yy(&mut self.yy, &self.state, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.boundaries
            .assemble_yy(&mut self.yy, &self.state, &self.ignored_eqs)?;
        Ok(())
    }

    pub fn calc_ff(&mut self) -> Result<(), StrError> {
        // clear vector
        self.ff.fill(0.0);

        // calculate all element local vectors
        let t = self.state.time;
        self.elements.assemble_ff(&mut self.ff, t, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.boundaries.assemble_ff(&mut self.ff, t, &self.ignored_eqs)?;

        // add concentrated loads
        for (eq, f) in &self.conc_loads {
            self.ff[*eq] += (f)(t);
        }
        Ok(())
    }

    /// Calculates Y (internal forces)
    pub fn calc_yy_to_delete(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // clear vector
        self.ls.yy.fill(0.0);

        // calculate all element local vectors
        self.elements.assemble_yy(&mut self.ls.yy, state, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.boundaries.assemble_yy(&mut self.ls.yy, state, &self.ignored_eqs)?;
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
    pub fn calc_ff_and_ddff_to_delete(&mut self, time: f64) -> Result<bool, StrError> {
        // make a copy of F and ΔF
        vec_copy(&mut self.ls.ff_old, &self.ls.ff).unwrap();
        vec_copy(&mut self.ls.ddff_old, &self.ls.ddff).unwrap();

        // update F ---------------------------------------------------------------

        // clear vector
        self.ls.ff.fill(0.0);

        // calculate all element local vectors
        self.elements.assemble_ff(&mut self.ls.ff, time, &self.ignored_eqs)?;

        // calculate all boundary elements local vectors
        self.boundaries.assemble_ff(&mut self.ls.ff, time, &self.ignored_eqs)?;

        // add concentrated loads
        for (eq, f) in &self.conc_loads {
            self.ls.ff[*eq] += (f)(time);
        }

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
        self.elements
            .assemble_kk_to_delete(&mut self.ls.kk, state, &self.ignored_eqs)?;
        self.boundaries
            .assemble_kk_to_delete(&mut self.ls.kk, state, &self.ignored_eqs)?;
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
            let val = self.presc_values[ip](state.time);
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
            Some("cannot start simulation because config.validate() failed")
        );
    }
}
