use super::{ElementsBoundary, ElementsInterior, FemState, OutputFiles};
use crate::base::{BcEssential, BcNatural, Config, Dof, Schema};
use crate::material::LocalState;
use crate::StrError;
use gemlab::mesh::{CellId, Mesh, PointId};
use russell_lab::{Stopwatch, Vector};
use russell_pde::EquationHandler;
use russell_sparse::{CooMatrix, Sym};
use std::collections::HashMap;
use std::fmt::Write;
use std::sync::Arc;
use uuid::Uuid;

/// Implements common (shared) functionality for all FEM solvers
pub struct FemData<'a> {
    /// Holds a unique identifier for this instance such that it can be tracked externally
    pub(crate) uuid: Uuid,

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

    pub(crate) state: FemState,
    pub(crate) ndof: usize,
    pub(crate) nu: usize,
    pub(crate) np: usize,
    pub(crate) nsys: usize,
    pub(crate) sym: Sym,
    pub(crate) nnz_kk: usize,
    pub(crate) nnz_kk_bar: usize,
    pub(crate) kk_check: CooMatrix,

    /// Pᵤ(t), lambda-free part of the prescribed values Ǔ = λ Pᵤ(t)
    pub(crate) ppu: Vector,
}

impl<'a> FemData<'a> {
    /// Allocates a new instance
    pub fn new(
        mesh: &Mesh,
        schema: &'a Schema,
        config: &'a Config,
        ebc: &'a BcEssential,
        nbc: &'a BcNatural,
    ) -> Result<Self, StrError> {
        // Check
        if let Some(msg) = config.validate() {
            println!("ERROR: {}", msg);
            return Err("cannot start simulation because config.validate() failed");
        }

        // Start stopwatch
        let stopwatch = Stopwatch::new();

        // Generate the list of prescribed equations and a map from equation to (PointId, Dof)
        let n_prescribed = ebc.functions.len();
        let mut p_list = Vec::with_capacity(n_prescribed);
        let mut eq_to_dof = HashMap::with_capacity(n_prescribed);
        for (point_id, dof) in ebc.functions.keys() {
            let eq = schema.dof_number(*point_id, *dof)?;
            p_list.push(eq);
            eq_to_dof.insert(eq, (*point_id, *dof));
        }

        // Allocate the equations handler
        let ndof = schema.ndof()?;
        let mut eq_handler = EquationHandler::new(ndof);
        eq_handler.recompute(&p_list);

        // Allocate array of functions to calculate prescribed values
        let mut presc_values = Vec::with_capacity(n_prescribed);
        for eq in eq_handler.prescribed() {
            let point_dof = eq_to_dof.get(eq).unwrap();
            let f = ebc.functions.get(point_dof).unwrap();
            presc_values.push(f.clone());
        }

        // Allocate array of concentrated loads
        let mut conc_loads = Vec::with_capacity(nbc.at_points.len());
        for (point_id, pbc, f) in &nbc.at_points {
            let eq = schema.dof_number(*point_id, pbc.dof())?;
            conc_loads.push((eq, f.clone()));
        }

        // Allocate auxiliary instances
        let boundaries = ElementsBoundary::new(mesh, schema, config, nbc)?;
        let mut elements = ElementsInterior::new(mesh, schema, config)?;

        // Array to ignore prescribed equations when building the reduced system
        let neq = eq_handler.neq(); // number of DOFs (without Lagrange multipliers)
        let mut ignored_eqs = vec![false; neq];
        if !config.lagrange_mult_method {
            for eq in eq_handler.prescribed() {
                ignored_eqs[*eq] = true;
            }
        };

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        let mut state = FemState::new(&mesh, &schema, &config)?;

        // Initialize internal variables
        elements.initialize_internal_values(&mut state)?;

        // Determine if the global stiffness matrix is symmetric and it's enabled
        let symmetric = if config.ignore_symmetry {
            false
        } else {
            elements.all_sym_kk() && boundaries.all_sym_kk()
        };

        // Determine symmetry type of the global stiffness matrix
        let genie = config.lin_sol_genie;
        let sym = genie.get_sym(symmetric);

        // Determine the system dimension
        let nu = eq_handler.nu();
        let np = eq_handler.np();
        let nsys = if config.lagrange_mult_method { ndof + np } else { nu };

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

        let yy = Vector::new(ndof);
        let ff = Vector::new(ndof);

        //////////////////////////////////////////////////////////////////////////////////////////////////////

        // Allocate output files handler
        let mut files = OutputFiles::new(mesh, schema, config, np)?;

        // Perform the first output
        files.execute(&schema, &config, &state, &yy)?;

        // return new instance
        Ok(FemData {
            uuid: Uuid::new_v4(),
            schema,
            config,
            eq_handler,
            presc_values,
            conc_loads,
            boundaries,
            elements,
            files,
            stopwatch,
            yy,
            ff,
            //
            state,
            ndof,
            nu,
            np,
            nsys,
            sym,
            nnz_kk,
            nnz_kk_bar,
            kk_check,
            ppu: Vector::new(np),
        })
    }

    /// Returns the number of equations effectively considered in the nonlinear system
    ///
    /// This number may be greater than the number of equations when using the Lagrange Multiplier Method,
    /// or this number may be smaller than the number of equations when using the System Partitioning Strategy.
    pub fn get_nsys(&self) -> usize {
        self.nsys
    }

    /// Returns the total number of degrees of freedom
    ///
    /// This number corresponds to the total number of equations without Lagrange multipliers and without prescribed values.
    pub fn get_ndof(&self) -> usize {
        self.ndof
    }

    /// Returns the index of a u component in the system used by the nonlinear solver (this is not U)
    pub fn get_u_index(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        let eq = self.schema.dof_number(point_id, dof)?;
        if self.config.lagrange_mult_method {
            Ok(eq)
        } else {
            Ok(self.eq_handler.iu(eq))
        }
    }

    pub fn get_state(&self) -> &FemState {
        &self.state
    }

    pub fn reset_algorithmic_variables(&mut self, load_reversal: bool) {
        self.state.reverse = load_reversal;
        self.elements.reset_algorithmic_variables(&mut self.state);
        self.state.reverse = false;
    }

    /// Returns the real simulation times (time) or loading increments (lambda)
    ///
    /// Time is used for transient/dynamic analyses, while lambda is used for steady/static analyses
    pub fn stations(&self) -> &Vec<f64> {
        self.files.stations()
    }

    /// Returns the history (time or lambda) of U components at selected points
    pub fn history_uu_comp(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        self.files.history_uu_comp(point_id, dof, &self.schema)
    }

    /// Returns the history (time or lambda) of Y (internal forces) components at selected points
    pub fn history_yy_comp(&self, point_id: PointId, dof: Dof) -> Option<&Vec<f64>> {
        self.files.history_yy_comp(point_id, dof, &self.schema)
    }

    /// Returns the history (time or lambda) of flux vectors at selected integration points
    pub fn history_local_fluxes(&self, cell_id: CellId) -> Option<&Vec<Vector>> {
        self.files.history_local_flux(cell_id)
    }

    /// Returns the history (time or lambda) of LocalState at selected integration points
    pub fn history_local_state(&self, cell_id: CellId) -> Option<&Vec<LocalState>> {
        self.files.history_local_state(cell_id)
    }

    /// Prints information about the system
    pub fn print_system_info(&self, continuation: &str) {
        if self.config.verbose {
            let mut b = vec![vec![String::new(); 3]; 3];
            let handler = if self.config.lagrange_mult_method { "LMM" } else { "SPS" };
            write!(&mut b[0][0], "neq  = {:?}", self.ndof).unwrap();
            write!(&mut b[1][0], "np   = {:?}", self.np).unwrap();
            write!(&mut b[2][0], "ndim = {:?}", self.nsys).unwrap();
            write!(&mut b[0][1], "nnz(K)     = {:?}", self.nnz_kk).unwrap();
            write!(&mut b[1][1], "nnz(K-bar) = {:?}", self.nnz_kk_bar).unwrap();
            write!(&mut b[2][1], "sym(K)     = {:?}", self.sym).unwrap();
            write!(&mut b[0][2], "genie        = {:?}", self.config.lin_sol_genie).unwrap();
            write!(&mut b[1][2], "continuation = {}", continuation).unwrap();
            write!(&mut b[2][2], "EBC handler  = {}", handler).unwrap();
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

    /// Initializes the nonlinear solver unknowns vector `u` from the state's `U`
    pub(crate) fn initialize_u(&self, u: &mut Vector) {
        if self.config.lagrange_mult_method {
            for eq in 0..self.ndof {
                u[eq] = self.state.uu[eq];
            }
        } else {
            for iu in 0..self.nu {
                let eq = self.eq_handler.unknown()[iu];
                u[iu] = self.state.uu[eq];
            }
        }
    }

    /// Calculates Pᵤ(t), lambda-free part of the prescribed values Ǔ = λ Pᵤ(t)
    pub(crate) fn calc_ppu(&mut self) {
        for ip in 0..self.np {
            self.ppu[ip] = self.presc_values[ip](self.state.time);
        }
    }

    /// Sets the state given the nonlinear solver variables (λ, u)
    ///
    /// This function requires that Cᵤ(t) (prescribed values) be calculated already.
    pub(crate) fn set_state(&mut self, l: f64, u: &Vector) {
        self.state.lambda = l;
        if self.config.lagrange_mult_method {
            for i in 0..self.ndof {
                self.state.uu[i] = u[i];
            }
        } else {
            for iu in 0..self.nu {
                let eq = self.eq_handler.unknown()[iu];
                self.state.uu[eq] = u[iu];
            }
            for ip in 0..self.np {
                let eq = self.eq_handler.prescribed()[ip];
                self.state.uu[eq] = l * self.ppu[ip];
            }
        }
    }

    /// Calculates Y, internal forces
    pub(crate) fn calc_yy(&mut self) -> Result<(), StrError> {
        // clear vector
        self.yy.fill(0.0);

        // calculate all element local vectors
        self.elements.assemble_yy(&mut self.yy, &self.state)?;

        // calculate all boundary elements local vectors
        self.boundaries.assemble_yy(&mut self.yy, &self.state)?;
        Ok(())
    }

    /// Calculates F(t), external forces
    pub(crate) fn calc_ff(&mut self) -> Result<(), StrError> {
        // clear vector
        self.ff.fill(0.0);

        // calculate all element local vectors
        let t = self.state.time;
        self.elements.assemble_ff(&mut self.ff, t)?;

        // calculate all boundary elements local vectors
        self.boundaries.assemble_ff(&mut self.ff, t)?;

        // add concentrated loads
        for (eq, f) in &self.conc_loads {
            self.ff[*eq] += (f)(t);
        }
        Ok(())
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
        let ebc = BcEssential::new();
        let nbc = BcNatural::new();

        // error due to config.validate
        let mut config = Config::new(&mesh);
        config.set_theta(0.0);
        assert_eq!(
            FemData::new(&mesh, &schema, &config, &ebc, &nbc).err(),
            Some("cannot start simulation because config.validate() failed")
        );
    }
}
