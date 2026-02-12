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

/// Holds the main data structures for the FEM simulation
pub struct FemData<'a> {
    /// Holds a unique identifier for this instance such that it can be tracked externally
    pub(crate) uuid: Uuid,

    /// Holds element types, material parameters, and specifies the DOF numbering schema
    pub(crate) schema: &'a Schema,

    /// Holds the configuration
    pub(crate) config: &'a Config<'a>,

    /// Stopwatch to measure computer time
    pub(crate) stopwatch: Stopwatch,

    /// Manages equation numbers (prescribed versus unknown)
    pub(crate) eq_handler: EquationHandler,

    /// Holds the functions to calculate the prescribed values
    ///
    /// Use `eq_handler.ip()` to access an entry in this array
    ///
    /// (np)
    pub(crate) presc_values: Vec<Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>>,

    /// Holds pairs of (dof_num, fn) to calculate concentrated loads
    pub(crate) conc_loads: Vec<(usize, Arc<dyn Fn(f64) -> f64 + Send + Sync + 'a>)>,

    // Holds a collection of boundary elements
    pub(crate) boundaries: ElementsBoundary<'a>,

    /// Holds a collection of elements
    pub(crate) elements: ElementsInterior<'a>,

    /// Number of degrees of freedom
    pub(crate) ndof: usize,

    /// Number of unknowns
    pub(crate) nu: usize,

    /// Number of prescribed DOFs
    pub(crate) np: usize,

    /// Number of equations in the nonlinear system (system dimension)
    pub(crate) nsys: usize,

    /// Symmetry type of the global stiffness matrix
    pub(crate) sym: Sym,

    /// Number of non-zero entries in the augmented LMM matrix M
    pub(crate) nnz_mm: usize,

    /// Number of non-zero entries in the SPS matrix K-bar
    pub(crate) nnz_kk_bar: usize,

    /// K-check matrix for the SPS method
    pub(crate) kk_check: CooMatrix,

    /// Holds the current state of the simulation
    pub(crate) state: FemState,

    /// Vector of internal forces
    ///
    /// (ndof)
    pub(crate) yy: Vector,

    /// Vector of external forces
    ///
    /// (ndof)
    pub(crate) ff: Vector,

    /// Pᵤ(t), lambda-free part of the prescribed values Ǔ = λ Pᵤ(t)
    ///
    /// (np)
    pub(crate) ppu: Vector,

    /// Handles output files
    pub(crate) files: OutputFiles,
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

        // Generate the list of prescribed DOFs and a map from DOF number to (PointId, Dof)
        let np = ebc.functions.len();
        let mut p_list = Vec::with_capacity(np);
        let mut i_to_dof = HashMap::with_capacity(np);
        for (point_id, dof) in ebc.functions.keys() {
            let i = schema.dof_number(*point_id, *dof)?;
            p_list.push(i);
            i_to_dof.insert(i, (*point_id, *dof));
        }

        // Allocate the equations handler
        let ndof = schema.ndof()?;
        let mut eq_handler = EquationHandler::new(ndof);
        eq_handler.recompute(&p_list);

        // Allocate array of functions to calculate prescribed values
        let mut presc_values = Vec::with_capacity(np);
        for i in eq_handler.prescribed() {
            let point_dof = i_to_dof.get(i).unwrap();
            let f = ebc.functions.get(point_dof).unwrap();
            presc_values.push(f.clone());
        }

        // Allocate array of concentrated loads
        let mut conc_loads = Vec::with_capacity(nbc.at_points.len());
        for (point_id, pbc, f) in &nbc.at_points {
            let i = schema.dof_number(*point_id, pbc.dof())?;
            conc_loads.push((i, f.clone()));
        }

        // Allocate elements
        let mut boundaries = ElementsBoundary::new(mesh, schema, config, nbc)?;
        let mut elements = ElementsInterior::new(mesh, schema, config)?;

        // Determine if the global stiffness matrix is symmetric
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
        let nsys = if config.lagrange_mult_method { ndof + np } else { nu };

        // Calculate the number of non-zero entries in the global stiffness matrix
        let mut nnz_mm = 0;
        let mut nnz_kk_bar = 0;
        let mut nnz_kk_check = 0;
        if config.lagrange_mult_method {
            elements.add_nnz_lmm(&mut nnz_mm, sym);
            boundaries.add_nnz_lmm(&mut nnz_mm, sym);
            if sym.triangular() {
                nnz_mm += np;
            } else {
                nnz_mm += 2 * np;
            }
        } else {
            elements.add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &eq_handler);
            boundaries.add_nnz_sps(&mut nnz_kk_bar, &mut nnz_kk_check, sym, &eq_handler);
        }

        // Allocate K-check matrix for SPS
        let kk_check = if config.lagrange_mult_method || np == 0 {
            CooMatrix::new(1, 1, 1, Sym::No).unwrap() // empty
        } else {
            CooMatrix::new(nu, np, nnz_kk_check, Sym::No).unwrap()
        };

        // Allocate the state
        let mut state = FemState::new(&mesh, &schema, &config)?;

        // Initialize internal variables
        elements.initialize_internal_values(&mut state)?;

        // Allocate Y, F, and Pᵤ vectors
        let mut yy = Vector::new(ndof);
        let ff = Vector::new(ndof);
        let ppu = Vector::new(np);

        // Calculate the first Y (only needed for output files)
        if config.out_history_yy_comp.len() > 0 {
            elements.assemble_yy(&mut yy, &state)?;
            boundaries.assemble_yy(&mut yy, &state)?;
        }

        // Allocate output files handler
        let mut files = OutputFiles::new(mesh, schema, config, np)?;

        // Perform the first output
        files.execute(&schema, &config, &state, &yy)?;

        // return new instance
        Ok(FemData {
            uuid: Uuid::new_v4(),
            schema,
            config,
            stopwatch,
            eq_handler,
            presc_values,
            conc_loads,
            boundaries,
            elements,
            ndof,
            nu,
            np,
            nsys,
            sym,
            nnz_mm,
            nnz_kk_bar,
            kk_check,
            state,
            yy,
            ff,
            ppu,
            files,
        })
    }

    /// Returns the total number of degrees of freedom (DOF)
    pub fn ndof(&self) -> usize {
        self.ndof
    }

    /// Returns the number of equations in the nonlinear system (system dimension)
    ///
    /// This number may be greater than the number of DOFs when using the Lagrange Multiplier Method,
    /// or it may be smaller than the number of DOFs when using the System Partitioning Strategy.
    pub fn nsys(&self) -> usize {
        self.nsys
    }

    /// Returns the index of an equation in the nonlinear system given a (PointId, Dof) pair
    ///
    /// If using the Lagrange Multiplier Method, the returned index is the same as the DOF number.
    /// If using the System Partitioning Strategy, the returned index corresponds to the index of an unknown.
    ///
    /// An error is returned if the (PointId, Dof) pair does not correspond to an unknown DOF.
    pub fn sys_index(&self, point_id: PointId, dof: Dof) -> Result<usize, StrError> {
        let i = self.schema.dof_number(point_id, dof)?;
        if self.config.lagrange_mult_method {
            Ok(i)
        } else {
            if self.eq_handler.is_unknown(i) {
                Ok(self.eq_handler.iu(i))
            } else {
                Err("the specified (PointId, Dof) pair does not correspond to an unknown DOF")
            }
        }
    }

    /// Returns an access the current state
    pub fn state(&self) -> &FemState {
        &self.state
    }

    /// Resets the algorithmic variables of all elements
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
            let genie = format!("{:?}", self.config.lin_sol_genie).to_ascii_uppercase();
            write!(&mut b[0][0], "ndof = {:?}", self.ndof).unwrap();
            write!(&mut b[1][0], "np   = {:?}", self.np).unwrap();
            write!(&mut b[2][0], "nsys = {:?}", self.nsys).unwrap();
            write!(&mut b[0][1], "nnz(M)     = {:?}", self.nnz_mm).unwrap();
            write!(&mut b[1][1], "nnz(K-bar) = {:?}", self.nnz_kk_bar).unwrap();
            write!(&mut b[2][1], "symmetry   = {:?}", self.sym).unwrap();
            write!(&mut b[0][2], "genie        = {}", genie).unwrap();
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
    pub(crate) fn initialize_sys_u(&self, u: &mut Vector) -> Result<(), StrError> {
        if self.config.lagrange_mult_method {
            for eq in 0..self.ndof {
                u[eq] = self.state.uu[eq];
            }
            if let Some(mu) = self.state.lag_mult.as_ref() {
                if mu.dim() != self.np {
                    return Err("The recorded Lagrange multipliers vector must have dimension equal to the actual number of prescribed values");
                }
                for ip in 0..self.np {
                    let j = self.ndof + ip;
                    u[j] = mu[ip];
                }
            }
        } else {
            for iu in 0..self.nu {
                let eq = self.eq_handler.unknown()[iu];
                u[iu] = self.state.uu[eq];
            }
        }
        Ok(())
    }

    // Records the Lagrange multipliers for future simulations
    pub(crate) fn record_lagrange_multipliers(&mut self, u: &Vector) {
        if self.config.lagrange_mult_method {
            let mut mu = Vector::new(self.np);
            for ip in 0..self.np {
                let j = self.ndof + ip;
                mu[ip] = u[j];
            }
            self.state.lag_mult = Some(mu);
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
        config.theta(0.0);
        assert_eq!(
            FemData::new(&mesh, &schema, &config, &ebc, &nbc).err(),
            Some("cannot start simulation because config.validate() failed")
        );
    }
}
