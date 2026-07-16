use super::{callback_history_e, callback_history_ep, callback_intersect, callback_ode_e, callback_ode_ep};
use super::{ep_jacobian, ep_residual};
use super::{ArgsExp, ArgsImp, LocalState, PlasticityTrait, PlotterData, Settings, StressStrainTrait};
use super::{CHEBYSHEV_TOL, F_TOL, HISTORY_N_OUT, PSEUDO_TIME_TOL};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{mat_inverse, mat_vec_mul};
use russell_lab::{InterpChebyshev, Matrix, NewtonSolver, RootFinder, Vector};
use russell_ode::{OdeSolver, Output, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, t4_ddot_t2_update};
use russell_tensor::{Tensor2, Tensor4};

/// Indicates the yield surface crossing case
#[derive(Clone, Copy, Debug)]
enum Case {
    AE,       // elastic
    AXB(f64), // elastic-elastoplastic; holds t_intersection
    BE,       // elastic; going inside (with eventual crossing)
    BXP(f64), // elastic-elastoplastic; going inside then outside with two crossings; holds t_intersection
    BP,       // elastoplastic
}

//// Holds the data for the explicit stress update algorithm
struct DataExp<'a> {
    /// Holds the arguments for the explicit stress update algorithm
    args: ArgsExp,

    /// Holds the solver for finding the yield surface intersection
    ode_intersection: OdeSolver<'a, ArgsExp>,

    /// Holds the solver for the elastic update
    ode_elastic: OdeSolver<'a, ArgsExp>,

    /// Holds the solver for the elastoplastic update
    ode_elastoplastic: OdeSolver<'a, ArgsExp>,

    /// Holds the ODE vector of unknowns for elastic case
    ode_y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    ode_y_ep: Vector,

    /// Holds the output during the intersection finding
    out_intersection: Output<'a, ArgsExp>,

    /// Holds the output during the elastic path
    out_history_el: Output<'a, ArgsExp>,

    /// Holds the output during the elastoplastic path
    out_history_ep: Output<'a, ArgsExp>,

    /// Holds the interpolant for finding the yield surface intersection
    interpolant: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// Enables recording stress-strain history
    save_history: bool,

    /// Holds the last Case analyzed by update_stress (for debugging)
    last_case: Option<Case>,
}

impl<'a> DataExp<'a> {
    /// Allocate a new instance
    fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // Allocate the interpolant for the explicit update
        let interp_nn_max = settings.gp_interp_nn_max();
        let interpolant = InterpChebyshev::new(interp_nn_max, 0.0, 1.0).unwrap();

        // Allocate interior stations for dense output (intersection finding)
        let chebyshev_points = InterpChebyshev::points(interp_nn_max);
        let interp_npoint = chebyshev_points.dim();
        let mut interior_t_out = vec![0.0; interp_npoint - 2];
        let xx_interior = &chebyshev_points.as_data()[1..(interp_npoint - 1)];
        xx_interior.into_iter().enumerate().for_each(|(i, x)| {
            interior_t_out[i] = (1.0 + x) / 2.0;
        });

        // Allocate arguments for the callback functions
        let args = ArgsExp::new(ideal, param, settings, interp_npoint)?;

        // ODE system: dσ/dt = Dₑ : Δε
        let ode_system_e = System::new(args.ndim_e, callback_ode_e);

        // ODE system: dσ/dt = Dₑₚ : Δε and dz/dt = λ h(σ,z)
        let ode_system_ep = System::new(args.ndim_ep, callback_ode_ep);

        // ODE solvers
        let ode_param = Params::new(settings.gp_ode_method());
        let ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        // Set function to handle yield surface intersection
        let mut out_intersection = Output::new();
        out_intersection
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(callback_intersect);

        // Set function to record the stress-strain history
        let mut out_history_el = Output::new();
        let mut out_history_ep = Output::new();
        let save_history = settings.gp_save_history();
        if save_history {
            let h_out = 1.0 / ((HISTORY_N_OUT - 1) as f64);
            out_history_el
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(callback_history_e);
            out_history_ep
                .set_dense_h_out(h_out)
                .unwrap()
                .set_dense_callback(callback_history_ep);
        }

        // Allocate ODE vectors
        let ode_y_e = Vector::new(args.ndim_e);
        let ode_y_ep = Vector::new(args.ndim_ep);

        // Allocate root finder
        let root_finder = RootFinder::new();

        // Done
        Ok(DataExp {
            args,
            ode_intersection,
            ode_elastic,
            ode_elastoplastic,
            ode_y_e,
            ode_y_ep,
            out_intersection,
            out_history_el,
            out_history_ep,
            interpolant,
            root_finder,
            save_history,
            last_case: None,
        })
    }
}

/// Holds the data for the implicit stress update algorithm
struct DataImp {
    /// Holds the arguments for the implicit stress update algorithm
    args: ArgsImp,

    /// Vector of unknowns for the local Newton-Raphson solver
    ///
    /// x := [σ, z, λ]
    x_newton: Vector,

    /// Jacobian matrix for the local Newton-Raphson solver
    jac_newton: Matrix,

    /// Inverse Jacobian matrix for the consistent tangent stiffness
    inv_jac_newton: Matrix,
}

impl<'a> DataImp {
    /// Allocate a new instance
    fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // Allocate arguments for the callback functions
        let args = ArgsImp::new(ideal, param, settings)?;

        // Allocate vector and matrix for the local Newton-Raphson solver (implicit integration)
        let ndim_nw = args.ncp + args.niv + 1; // dimension of the local nonlinear problem r = [rσ, rz, rλ] = 0
        let x_newton = Vector::new(ndim_nw);
        let jac_newton = Matrix::new(ndim_nw, ndim_nw);
        let inv_jac_newton = Matrix::new(ndim_nw, ndim_nw);

        // Done
        Ok(DataImp {
            args,
            x_newton,
            jac_newton,
            inv_jac_newton,
        })
    }
}

/// Implements general elastoplasticity models
pub struct Elastoplastic<'a> {
    /// Holds the data for the explicit stress update algorithm
    data_exp: Option<DataExp<'a>>,

    /// Holds the data for the implicit stress update algorithm
    data_imp: Option<DataImp>,

    /// Enables the explicit stress-update
    explicit_update: bool,

    /// Enables verbose mode
    verbose: bool,
}

impl<'a> Elastoplastic<'a> {
    /// Returns a reference to the model
    fn model_ref(&self) -> &dyn PlasticityTrait {
        if self.explicit_update {
            self.data_exp.as_ref().unwrap().args.model.as_ref()
        } else {
            self.data_imp.as_ref().unwrap().args.model.as_ref()
        }
    }
    
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        let (data_exp, data_imp) = if settings.gp_explicit_update() {
            (Some(DataExp::new(ideal, param, settings)?), None)
        } else {
            (None, Some(DataImp::new(ideal, param, settings)?))
        };
        // done
        Ok(Elastoplastic {
            data_exp,
            data_imp,
            explicit_update: settings.gp_explicit_update(),
            verbose: false,
        })
    }

    /// Calculates the yield function f
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        if let Some(data) = self.data_exp.as_ref() {
            data.args.model.calc_f(state)
        } else {
            let data = self.data_imp.as_ref().unwrap();
            data.args.model.calc_f(state)
        }
    }

    /// Returns the stress-strain history during the intersection finding (e.g., for debugging)
    pub fn get_history_int(&self) -> Result<PlotterData, StrError> {
        if !self.explicit_update {
            return Err("history is only available for explicit update");
        }
        match self.data_exp.as_ref().unwrap().args.history_int.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub fn get_history_eep(&self) -> Result<PlotterData, StrError> {
        if !self.explicit_update {
            return Err("history is only available for explicit update");
        }
        match self.data_exp.as_ref().unwrap().args.history_eep.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns true if the trial stress path leads to the inside of the yield surface
    ///
    /// Warning: this function must only be called if the stress point is on (or near) the yield surface.
    fn going_inside(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<bool, StrError> {
        // check that the explicit update is enabled
        if !self.explicit_update {
            return Err("going_inside is only available for explicit update");
        }

        // gradients of the yield function
        let data = self.data_exp.as_mut().unwrap();
        data.args.model.calc_fs(&mut data.args.fs, state)?;

        // Dₑ
        data.args.model.calc_dde(&mut data.args.dde, state)?;

        // (df/dσ) : Dₑ : Δε
        let indicator = t2_ddot_t4_ddot_t2(&data.args.fs, &data.args.dde, delta_strain);
        Ok(indicator < 0.0)
    }

    /// Performs the intersection finding algorithm
    ///
    /// Returns `(t_int, yf_trial)`
    fn intersection_finding(&mut self, state: &LocalState, inside: bool) -> Result<(Option<f64>, f64), StrError> {
        // check that the explicit update is enabled
        if !self.explicit_update {
            return Err("intersection_finding is only available for explicit update");
        }

        // copy z into arguments (z is frozen)
        let data = self.data_exp.as_mut().unwrap();
        data.args.state.int_vars.set_vector(state.int_vars.as_data());

        // copy σ into {y}
        data.ode_y_e.set_vector(state.stress.vector().as_data());

        // solve the elastic problem with intersection finding data
        data.ode_intersection.solve(
            &mut data.ode_y_e,
            0.0,
            1.0,
            None,
            &mut data.args,
            Some(&mut data.out_intersection),
        )?;
        assert_eq!(data.args.yf_count, data.args.yf_values.dim());

        // set data for interpolation
        data.interpolant
            .adapt_data(CHEBYSHEV_TOL, data.args.yf_values.as_data())?;

        // find roots == intersections
        let roots = data.root_finder.chebyshev(&data.interpolant)?;

        // extract last root (ignore first root if crossing twice)
        let t_int = if inside {
            match roots.len() {
                0 => None,
                1 => Some(roots[0]),
                _ => return Err("inside: cannot handle more than one intersection"),
            }
        } else {
            match roots.len() {
                0 => None,
                1 => None,
                2 => Some(roots[1]),
                _ => return Err("not inside: cannot handle more than two intersections"),
            }
        };

        // trial yield function value
        let yf_trial = data.args.yf_values[data.args.yf_count - 1];

        // results
        Ok((t_int, yf_trial))
    }

    /// Selects the yield surface crossing case
    fn select_case(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<Case, StrError> {
        // check that the explicit update is enabled
        if !self.explicit_update {
            return Err("select_case is only available for explicit update");
        }

        // current yield function value: f(σ, z)
        let yf_initial = self.data_exp.as_ref().unwrap().args.model.calc_f(state)?;

        // run analysis
        if yf_initial < 0.0 {
            //
            // A: inside the yield surface
            //
            let (t_intersection, yf_trial) = self.intersection_finding(state, true)?;
            match t_intersection {
                Some(t_int) => {
                    if t_int <= PSEUDO_TIME_TOL {
                        // start inside, intersecting the YS with a tiny length, meaning that
                        // the stress point is very close to the yield surface (from the inside)
                        // in this situation, disregard the elastic regime altogether => DH
                        Ok(Case::BP)
                    } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                        // start inside, intersecting the YS after crossing the "whole" elastic domain
                        Ok(Case::AE)
                    } else {
                        // start inside, crossing the yield surface
                        Ok(Case::AXB(t_int))
                    }
                }
                None => {
                    assert!(yf_trial <= 0.0); // cannot be positive if there is no intersection
                    Ok(Case::AE)
                }
            }
        } else {
            //
            // D: on the yield surface or slightly outside
            //
            if self.going_inside(state, delta_strain)? {
                let (t_intersection, yf_trial) = self.intersection_finding(state, false)?;
                match t_intersection {
                    Some(t_int) => {
                        if t_int <= PSEUDO_TIME_TOL {
                            // start on YS, going inside with a tiny length, meaning that
                            // the stress point remains on the yield surface due to a
                            // "neutral loading"
                            Ok(Case::BP)
                        } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                            // start on YS, going inside, reaching the "other" side of the YS
                            Ok(Case::BE)
                        } else {
                            // start on YS, crossing the "whole" elastic domain,
                            // and reaching the outside again
                            Ok(Case::BXP(t_int))
                        }
                    }
                    None => {
                        assert!(yf_trial <= 0.0); // cannot be positive if there is no intersection
                        Ok(Case::BE)
                    }
                }
            } else {
                Ok(Case::BP)
            }
        }
    }

    /// Calculates the consistent tangent stiffness for the explicit method (not available)
    fn explicit_stiffness(&mut self, _dd: &mut Tensor4, _state: &LocalState) -> Result<(), StrError> {
        Err("stiffness is not available for explicit update")
    }

    /// Updates the stress tensor given the strain increment tensor using the explicit method
    fn explicit_update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError> {
        {
            let data = self.data_exp.as_mut().unwrap();

            // set Δε in arguments struct
            data.args.del_eps.set_tensor(1.0, delta_strain);

            // enable history
            if data.save_history {
                match state.strain.as_ref() {
                    Some(strain) => data.args.state.strain = Some(strain.clone()),
                    None => {
                        return Err("state must have strain enabled");
                    }
                }
                data.args.history_int = Some(PlotterData::new());
                data.args.history_eep = Some(PlotterData::new());
            }
        }

        // select case regarding yield surface crossing
        let case = self.select_case(state, delta_strain)?;

        // get mutable reference to data
        let data = self.data_exp.as_mut().unwrap();

        // perform the update
        match case {
            // purely elastic => done
            Case::AE | Case::BE => {
                // update (note that select_case already calculated this path)
                state.stress.vector_mut().set_vector(data.ode_y_e.as_data());
                state.elastic = true;
            }

            // elastic-elastoplastic with crossing
            Case::AXB(t_int) | Case::BXP(t_int) => {
                // copy σ into {y} (again; to start from scratch; because select_case modified ode_y_e)
                data.ode_y_e.set_vector(state.stress.vector().as_data());

                // solve the elastic problem (again) to update σ to the intersection point
                data.ode_elastic.solve(
                    &mut data.ode_y_e,
                    0.0,
                    t_int,
                    None,
                    &mut data.args,
                    Some(&mut data.out_history_el),
                )?;

                // set stress at intersection
                state.stress.vector_mut().set_vector(data.ode_y_e.as_data());

                // elastoplastic run: join σ and z into {y} (now z plays a role)
                data.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());

                // solve elastoplastic problem (starting from t_int)
                data.ode_elastoplastic.solve(
                    &mut data.ode_y_ep,
                    t_int,
                    1.0,
                    None,
                    &mut data.args,
                    Some(&mut data.out_history_ep),
                )?;

                // update: split {y} into σ and z
                data.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }

            // elastoplastic
            Case::BP => {
                // join σ and z into {y} (now z plays a role)
                data.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.int_vars.as_data());

                // solve elastoplastic problem
                data.ode_elastoplastic.solve(
                    &mut data.ode_y_ep,
                    0.0,
                    1.0,
                    None,
                    &mut data.args,
                    Some(&mut data.out_history_ep),
                )?;

                // update: split {y} into σ and z
                data.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.int_vars.as_mut_data());
                state.elastic = false;
            }
        }

        // print message
        if self.verbose {
            println!("👉 {:?}", case);
        }

        // record last_case for debugging
        data.last_case = Some(case);
        Ok(())
    }

    /// Calculates the consistent tangent stiffness for the implicit method
    fn implicit_stiffness(&mut self, dd: &mut Tensor4, state: &LocalState) -> Result<(), StrError> {
        let data = self.data_imp.as_mut().unwrap();

        // Calculate the elastic moduli if they have not been calculated yet
        if !data.args.elastic_moduli_calculated {
            // Calculate Dₑ
            data.args.model.calc_dde(&mut data.args.dde, state)?;

            // Calculate Cₑ
            mat_inverse(data.args.cce.matrix_mut(), data.args.dde.matrix())?;

            // Set flag
            data.args.elastic_moduli_calculated = true;
        }

        // Handle elastic case
        if state.elastic {
            dd.set_tensor(1.0, &data.args.dde); // D ← Dₑ
            return Ok(());
        }

        // --- Elastoplastic stiffness ---

        // Set some auxiliary constants
        let ns = data.args.ncp; // number of stress components
        let nz = data.args.niv; // number of internal variables
        let nsz = ns + nz; // index of λ in x

        // Build vector of unknowns x := [σ, z, λ]
        for i in 0..ns {
            data.x_newton[i] = state.stress.vector()[i];
        }
        for i in 0..nz {
            data.x_newton[ns + i] = state.int_vars[i];
        }
        data.x_newton[nsz] = state.lambda_alg;

        // Calculate the Jacobian matrix
        ep_jacobian(&mut data.jac_newton, &data.x_newton, &mut data.args)?;

        // Invert the Jacobian matrix to get the consistent tangent stiffness
        mat_inverse(&mut data.inv_jac_newton, &data.jac_newton)?;

        // Set the consistent tangent stiffness: D ← = inv(Jacobian)ₛₛ
        for i in 0..ns {
            for j in 0..ns {
                dd.matrix_mut().set(i, j, data.inv_jac_newton.get(i, j));
            }
        }
        Ok(())
    }

    /// Updates the stress tensor given the strain increment tensor using the implicit method
    fn implicit_update_stress(&mut self, state: &mut LocalState, delta_strain: &Tensor2) -> Result<(), StrError> {
        let data = self.data_imp.as_mut().unwrap();

        // Calculate the elastic moduli if they have not been calculated yet
        if !data.args.elastic_moduli_calculated {
            // Calculate Dₑ
            data.args.model.calc_dde(&mut data.args.dde, state)?;

            // Calculate Cₑ
            mat_inverse(data.args.cce.matrix_mut(), data.args.dde.matrix())?;

            // Set flag
            data.args.elastic_moduli_calculated = true;
        }

        // Reset data to elastic state
        state.elastic = true; // aka, unloading
        state.lambda_alg = 0.0; // algorithmic Lagrange multiplier

        // Note that, at this stage:
        // 1. σ_old = σ_current = state.stress
        // 2. z_old = z_current = state.int_vars

        // Trial update: σ_trial = σ_old + Dₑ : Δε thus σ += Dₑ : Δε
        t4_ddot_t2_update(&mut state.stress, 1.0, &data.args.dde, delta_strain, 1.0);

        // Trial yield function value: f(σ_trial, z_old)
        let f_trial = data.args.model.calc_f(state)?;

        // Exit on elastic update
        if f_trial < F_TOL * data.args.model.calc_f_ref() {
            // Elastic update: σ = σ_trial, z = z_old, λ_alg = 0.0
            return Ok(());
        }

        // --- Elastoplastic update ---

        // Calculate ε_trial = Cₑ : σ_trial
        mat_vec_mul(
            &mut data.args.eps_trial,
            1.0,
            data.args.cce.matrix(),
            state.stress.vector(),
        )?;

        // Set z_old in arguments struct
        data.args.z_old.set_vector(state.int_vars.as_data());

        // Set some auxiliary constants
        let ns = data.args.ncp; // number of stress components
        let nz = data.args.niv; // number of internal variables
        let nsz = ns + nz; // index of λ in x

        // Build vector of unknowns x := [σ, z, λ]
        for i in 0..ns {
            data.x_newton[i] = state.stress.vector()[i];
        }
        for i in 0..nz {
            data.x_newton[ns + i] = state.int_vars[i];
        }
        data.x_newton[nsz] = state.lambda_alg; // initial guess

        // Solve the nonlinear system of equations
        let ndim = ns + nz + 1; // dimension of the local nonlinear problem r = 0
        let mut newton = NewtonSolver::new(ndim)?;
        newton.solve(&mut data.x_newton, &mut data.args, ep_residual, ep_jacobian)?;

        // Copy the results back into the state
        for i in 0..ns {
            state.stress.vector_mut()[i] = data.x_newton[i];
        }
        for i in 0..nz {
            state.int_vars[i] = data.x_newton[ns + i];
        }
        state.lambda_alg = data.x_newton[nsz];

        // Set the elastic flag to false (elastoplastic update)
        state.elastic = false;

        // Done
        Ok(())
    }
}

impl<'a> StressStrainTrait for Elastoplastic<'a> {
    /// Indicates that the stiffness matrix is symmetric
    fn symmetric_stiffness(&self) -> bool {
        self.model_ref().symmetric_stiffness()
    }

    /// Returns the number of internal variables
    fn n_int_vars(&self) -> usize {
        self.model_ref().n_int_vars()
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.model_ref().initialize_int_vars(state)
    }

    /// Computes the consistent tangent stiffness
    fn stiffness(
        &mut self,
        dd: &mut Tensor4,
        state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        if self.explicit_update {
            self.explicit_stiffness(dd, state)
        } else {
            self.implicit_stiffness(dd, state)
        }
    }

    /// Updates the stress tensor given the strain increment tensor
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        if self.explicit_update {
            self.explicit_update_stress(state, delta_strain)
        } else {
            self.implicit_update_stress(state, delta_strain)
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Case, Elastoplastic};
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::{extract_von_mises_params, extract_von_mises_params_kg};
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, StressStrainTrait, VonMises};
    use plotpy::Text;
    use russell_lab::{approx_eq, math::PI};
    use russell_lab::{mat_approx_eq, vec_approx_eq};
    use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};
    use std::collections::HashMap;

    const VERBOSE: bool = true;
    const SAVE_FIGURE: bool = false;

    // Returns a vector of keys associated with the Case (for debugging)
    fn case_to_keys(case: &Case) -> Vec<&str> {
        match case {
            Case::AE => vec!["A", "E"],
            Case::AXB(..) => vec!["A", "X", "B"],
            Case::BE => vec!["B", "E"],
            Case::BXP(..) => vec!["B", "X", "P"],
            Case::BP => vec!["B", "P"],
        }
    }

    // Generates a initial state for the von Mises model
    fn gen_ini_state_von_mises(
        ideal: &Idealization,  // geometry idealization
        model: &Elastoplastic, // model
        p: f64,                // initial mean invariant
        q: f64,                // initial deviatoric invariant (only if yf_error is None)
        alpha: f64,            // initial octahedral angle related to the Lode invariant
    ) -> LocalState {
        let distance = p * SQRT_3;
        let radius = q * SQRT_2_BY_3;
        let n_int_vars = model.n_int_vars();
        let mut state = LocalState::new(ideal.mandel(), n_int_vars);
        state.stress = Tensor2::new_from_octahedral_alpha(distance, radius, alpha, ideal.two_dim).unwrap();
        model.initialize_int_vars(&mut state).unwrap();
        state.enable_strain(); // for plotting
        state
    }

    // Runs the stress-update with the von Mises model
    //
    // returns (deps_v, deps_d)
    fn update_with_von_mises(
        param: &StressStrain,      // parameters
        model: &mut Elastoplastic, // model
        state: &mut LocalState,    // the state to be updated
        p_el: f64,                 // next mean stress corresponding to a linear elastic path
        q_el: f64,                 // next deviatoric stress corresponding to a linear elastic path
        alpha_el: f64,             // next octahedral angle corresponding to a linear elastic path
    ) -> (f64, f64) {
        // calculate stress increment
        let distance = p_el * SQRT_3;
        let radius = q_el * SQRT_2_BY_3;
        let mandel = state.stress.mandel();
        let two_dim = mandel.two_dim();
        let stress_fin = Tensor2::new_from_octahedral_alpha(distance, radius, alpha_el, two_dim).unwrap();
        let mut dsigma = Tensor2::new(mandel);
        t2_add(&mut dsigma, 1.0, &stress_fin, -1.0, &state.stress); // Δσ = σ_fin - σ_ini

        // calculate strain increment
        let (young, poisson, _, _) = extract_von_mises_params(param);
        let elast = LinElasticity::new(young, poisson, two_dim, false);
        let mut cc = Tensor4::new(mandel);
        elast.calc_compliance(&mut cc).unwrap();
        let mut depsilon = Tensor2::new(mandel);
        t4_ddot_t2(&mut depsilon, 1.0, &cc, &dsigma); // Δε = C : Δσ

        // perform the update
        model.update_stress(state, &depsilon, 0, 0).unwrap(); // update stress
        state.strain.as_mut().unwrap().update(1.0, &depsilon); // update strain (for plotting)
        (depsilon.invariant_eps_v(), depsilon.invariant_eps_d())
    }

    // Returns Text for labels in plots
    fn get_text_label() -> Text {
        let mut text = Text::new();
        text.set_fontsize(12.0)
            .set_bbox(true)
            .set_bbox_style("round,pad=0.1")
            .set_bbox_facecolor("#fff8c1")
            .set_bbox_edgecolor("#7a7a7a")
            .set_align_horizontal("center")
            .set_align_vertical("center");
        text
    }

    // Plot the results (test type # a)
    fn do_plot_a(
        file_stem: &str,
        data: &HashMap<i32, Vec<LocalState>>,
        labels_oct: &[(&str, f64, f64)],
        labels_tyf: &[(&str, f64, f64)],
        oct_radius_max: Option<f64>,
        tyf_range: Option<(f64, f64)>,
    ) {
        let mut plotter = Plotter::new();
        plotter.set_layout_selected_2x2(Axis::Time, Axis::Yield);
        if let Some(r) = oct_radius_max {
            plotter.set_oct_radius_max(r);
        }
        for (lode, marker, size, void) in [(-1, "s", 10.0, true), (0, "o", 8.0, true), (1, ".", 8.0, false)] {
            let states = data.get(&lode).unwrap();
            let mut data = PlotterData::new();
            for i in 0..states.len() {
                let s = &states[i];
                let f = s.stress.invariant_q() - s.int_vars[0];
                let t = (i as f64) / 2.0;
                data.push(&s.stress, s.strain.as_ref(), Some(f), Some(t));
            }
            plotter
                .add_2x2(&data, false, |curve, _, _| {
                    curve
                        .set_marker_style(marker)
                        .set_marker_size(size)
                        .set_marker_void(void)
                        .set_label(&format!(" $\\ell = {}$", lode));
                })
                .unwrap();
            if lode == 0 {
                let p = states.len() - 1;
                let radius_0 = states[0].int_vars[0] * SQRT_2_BY_3;
                let radius_1 = states[p].int_vars[0] * SQRT_2_BY_3;
                plotter.set_oct_circle(radius_0, |_| {});
                plotter.set_oct_circle(radius_1, |canvas| {
                    canvas.set_line_style("-");
                });
            }
        }
        plotter.set_extra(Axis::OctX, Axis::OctY, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_oct {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
        });
        plotter.set_extra(Axis::Time, Axis::Yield, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_tyf {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
            if let Some((y_min, y_max)) = tyf_range {
                plot.set_yrange(y_min, y_max);
            }
        });
        plotter.save(&format!("/tmp/pmsim/material/{}.svg", file_stem)).unwrap();
    }

    // Plot the results (test type # b)
    fn do_plot_b(
        file_stem: &str,
        model: &Elastoplastic,
        states: &[LocalState],
        labels_oct: &[(&str, f64, f64)],
        labels_tyf: &[(&str, f64, f64)],
        oct_radius_max: Option<f64>,
        tyf_range: Option<(f64, f64)>,
    ) {
        let mut plotter = Plotter::new();
        plotter
            .set_tab_leg_ncol(2)
            .set_layout_selected_2x2(Axis::Time, Axis::Yield);
        if let Some(r) = oct_radius_max {
            plotter.set_oct_radius_max(r);
        }
        let history_int = model.get_history_int().unwrap();
        let history_eep = model.get_history_eep().unwrap();
        plotter
            .add_2x2(&history_int, false, |curve, _, _| {
                curve
                    .set_label("history(int)")
                    .set_line_color("gold")
                    .set_line_style("-");
            })
            .unwrap();
        plotter
            .add_2x2(&history_eep, false, |curve, _, _| {
                curve
                    .set_label("history(e-ep)")
                    .set_line_color("#7a7a7a")
                    .set_line_style("--")
                    .set_marker_style(".")
                    .set_marker_every(2);
            })
            .unwrap();
        let mut data = PlotterData::new();
        for i in 0..states.len() {
            let s = &states[i];
            let f = model.yield_function(s).unwrap();
            let t = i as f64;
            data.push(&s.stress, s.strain.as_ref(), Some(f), Some(t));
        }
        plotter
            .add_2x2(&data, false, |curve, _, _| {
                curve
                    .set_label("actual update")
                    .set_marker_style("s")
                    .set_marker_void(true);
            })
            .unwrap();
        let p = states.len() - 1;
        let radius_0 = states[0].int_vars[0] * SQRT_2_BY_3;
        let radius_1 = states[p].int_vars[0] * SQRT_2_BY_3;
        plotter.set_oct_circle(radius_0, |_| {});
        plotter.set_oct_circle(radius_1, |canvas| {
            canvas.set_line_style("-");
        });
        plotter.set_extra(Axis::OctX, Axis::OctY, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_oct {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
        });
        plotter.set_extra(Axis::Time, Axis::Yield, move |plot| {
            let mut text = get_text_label();
            for (label, x, y) in labels_tyf {
                text.draw(*x, *y, label);
            }
            plot.add(&text);
            if let Some((y_min, y_max)) = tyf_range {
                plot.set_yrange(y_min, y_max);
            }
        });
        plotter.save(&format!("/tmp/pmsim/material/{}.svg", file_stem)).unwrap();
    }

    #[test]
    fn update_stress_von_mises_1() {
        //
        // This tests runs 2D and 3D stress updates with three Lode angles
        // First, the yield surface is reached exactly with one step.
        // Second, an elastoplastic update is induced.
        //
        // Cases: AB, AC, and DH

        // parameters
        let param = StressStrain::sample_von_mises();
        let mut settings = Settings::new();
        settings.set_gp_explicit_update(true);
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 2.0);
        let (sig_m_1, sig_d_1) = (1.0, z_ini); // will reach the yield surface exactly
        let (sig_m_2, sig_d_2) = (2.0, 2.0 * z_ini); // to calc the next elastic trial increment

        // data for plotting
        let mut data_2d = HashMap::new(); // map: lode => states (2D only)

        // test
        for ndim in [2, 3] {
            for lode_int in [-1, 0, 1] {
                let lode = lode_int as f64;
                let alpha = PI / 2.0 - f64::acos(lode) / 3.0;
                let alpha_deg = alpha * 180.0 / PI;
                if VERBOSE {
                    println!("\nndim = {}, lode = {}, alpha = {}°", ndim, lode, alpha_deg);
                }

                // model
                let ideal = Idealization::new(ndim);
                let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
                model.verbose = VERBOSE;

                // initial state
                let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
                if ndim == 2 {
                    data_2d.insert(lode_int, vec![state.clone()]);
                }

                // Cases AB or AC: elastic update (to yield surface exactly)
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha);
                let sig_m_1 = state.stress.invariant_p();
                let sig_d_1 = state.stress.invariant_q();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }

                // check
                let correct_sig_m = sig_m_0 + kk * deps_v;
                let correct_sig_d = sig_d_0 + 3.0 * gg * deps_d;
                approx_eq(sig_m_1, correct_sig_m, 1e-14);
                approx_eq(sig_d_1, correct_sig_d, 1e-13);
                approx_eq(state.int_vars[0], z_ini, 1e-15);
                assert_eq!(state.elastic, true);
                let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, ["A", "E"]);

                // Case DH: elastoplastic update
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_2, sig_d_2, alpha);
                let sig_m_2 = state.stress.invariant_p();
                let sig_d_2 = state.stress.invariant_q();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }

                // check
                let correct_sig_m = sig_m_1 + kk * deps_v;
                let correct_sig_d = sig_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
                approx_eq(sig_m_2, correct_sig_m, 1e-14);
                approx_eq(sig_d_2, correct_sig_d, 1e-13);
                approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
                assert_eq!(state.elastic, false);
                let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, &["B", "P"]);
            }
        }

        // plot
        if SAVE_FIGURE {
            let labels_oct = [
                ("A", 0.0, -2.3),
                ("E,B", 6.5, 1.5),
                ("E", -1.5, 7.1),
                ("E", 2.0, 6.5),
                ("P", 10.5, 4.0),
                ("P", 3.2, 9.5),
                ("P", -1.5, 10.0),
            ];
            let labels_tyf = [("A", 0.0, -7.9), ("E", 0.5, 1.1), ("B", 0.5, -1.1), ("P", 1.0, -1.1)];
            do_plot_a(
                "test_update_stress_von_mises_1",
                &data_2d,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_2() {
        //
        // This test simulates an elastoplastic update with the
        // stress point crossing the yield surface
        //
        // Case AXD

        // parameters
        let param = StressStrain::sample_von_mises();
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, z_ini + 9.0, PI / 3.0); // will cross the yield surface

        // settings
        let mut settings = Settings::new();
        settings.set_gp_explicit_update(true).set_gp_save_history(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());

        // check
        let deps_d_e = z_ini / (3.0 * gg);
        let deps_d_ep = deps_d - deps_d_e;
        let correct_sig_m = kk * deps_v;
        let correct_sig_d = z_ini + 3.0 * gg * hh * deps_d_ep / (3.0 * gg + hh);
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-13);
        approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
        assert_eq!(state.elastic, false);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["A", "X", "B"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("A", 0.0, -2.0), ("X", 4.9, 5.5), ("B", 6.8, 8.5)];
            let labels_tyf = [("A", 0.0, -8.0), ("X", 0.46, 0.71), ("B", 1.0, 0.9)];
            do_plot_b(
                "test_update_stress_von_mises_2",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(9.5),
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3a() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "other" side.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (0.0, 0.99999999999999999);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0);

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", -5.2, -7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3a",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3b() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "other" side.
        // Nonetheless, now an initial drift is induced making the stress update
        // algorithm to ignore the first yield surface crossing.
        //
        // Case DE

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1.0, 0.8);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0); // going inside, after crossing because of drift

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.7, 7.2), ("E", -1.0, -5.0)];
            let labels_tyf = [("B", 0.0, -0.12), ("E", 1.0, -0.8)];
            do_plot_b(
                "test_update_stress_von_mises_3b",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3c() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "lower" side,
        // according to a predefined "vertical" path on the octahedral plane.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (radius_0 * f64::cos(alpha_0), -radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", 5.2, -7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3c",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_3d() {
        //
        // This test simulates a purely elastic update with the stress point
        // starting from the yield surface and crossing to the "left" side,
        // according to a predefined "horizontal" path on the octahedral plane.
        //
        // Case DF

        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (-radius_0 * f64::cos(alpha_0), radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());

        // check
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.int_vars[0], z_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "E"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.2, 7.2), ("E", -5.2, 7.2)];
            let labels_tyf = [("B", 0.0, 1.0), ("E", 1.0, 1.0)];
            do_plot_b(
                "test_update_stress_von_mises_3d",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                None,
                Some((-10.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_4() {
        //
        // This test simulates an elastic-elastoplastic update with an initial drift.
        // The elastoplastic update happens after an initial elastic update before
        // the yield surface intersection is found.
        //
        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1.0, 2.5);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -PI / 3.0);

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());

        // check
        assert_eq!(state.elastic, false);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "X", "P"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 5.7, 7.2), ("X", 8.0, -3.0), ("P", 2.8, -10.0)];
            let labels_tyf = [("B", 0.0, 1.0), ("X", 0.41, 0.0), ("P", 1.0, 0.4)];
            do_plot_b(
                "test_update_stress_von_mises_4",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(10.5),
                Some((-3.0, 2.0)),
            );
        }
    }

    #[test]
    fn update_stress_von_mises_5() {
        //
        // This test simulates an elastoplastic loading such that the initial
        // loading direction is tangent to the yield surface. The initial drift
        // is essential to make the Case DH with going_inside = true to activate.
        // In this situation:
        //     yf_initial = 1.2434497875801753E-14
        //     indicator = -1.2434497875801753E-14 (going inside; but actually tangent)
        //     roots = [0, 0] (double root at the initial state)
        //     has_intersection = false
        //     yf_trial = 9
        //
        // parameters
        let param = StressStrain::sample_von_mises();
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);

        // constants
        let (drift, mz) = (1e-14, 2.0);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, 0.0);

        // settings
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);

        // model
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = Elastoplastic::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;

        // initial state
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);

        // array of states for plotting
        let mut states = vec![state.clone()];

        // update, crossing the yield surface
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());

        // check
        assert_eq!(state.elastic, false);
        let case = model.data_exp.as_ref().unwrap().last_case.as_ref().unwrap();
        let keys = case_to_keys(case);
        assert_eq!(keys, &["B", "P"]);

        // plot
        if SAVE_FIGURE {
            let labels_oct = [("B", 3.0, 8.5), ("P", 11.5, -2.0)];
            let labels_tyf = [("B", 0.0, -0.3), ("P", 1.0, -0.3)];
            do_plot_b(
                "test_update_stress_von_mises_5",
                &model,
                &states,
                &labels_oct,
                &labels_tyf,
                Some(9.5),
                Some((-2.0, 2.0)),
            );
        }
    }

    const YOUNG: f64 = 1500.0;
    const POISSON: f64 = 0.45;
    const HH: f64 = 800.0;
    const Z_INI: f64 = 9.0;

    #[test]
    fn implicit_consistent_modulus_matches_von_mises() {
        // Allocate von Mises model directly
        let ideal = Idealization::new(2);
        let param = StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            hh: HH,
            z_ini: Z_INI,
        };
        let settings = Settings::new();
        let mut vm = VonMises::new(&ideal, &param, &settings).unwrap();

        // Set the initial state (zero stress; inside yield surface)
        let mandel = ideal.mandel();
        let n_int_vars = vm.n_int_vars();
        let mut state0 = LocalState::new(mandel, n_int_vars);
        state0.enable_strain();
        vm.initialize_int_vars(&mut state0).unwrap();
        assert_eq!(state0.int_vars[0], Z_INI);

        // Allocate the Elastoplastic wrapper to von Mises model
        let mut ep = Elastoplastic::new(&ideal, &param, &settings).unwrap();

        // Calculate the consistent tangent modulus at the initial state using the von Mises model
        let mut dd_vm = Tensor4::new(mandel);
        vm.stiffness(&mut dd_vm, &state0, 0, 0).unwrap();
        // println!("dd_vm =\n{}", dd_vm.as_matrix());

        // Calculate the consistent tangent modulus at the initial state using the Elastoplastic wrapper
        let mut dd_ep = Tensor4::new(mandel);
        ep.stiffness(&mut dd_ep, &state0, 0, 0).unwrap();
        // println!("dd_ep =\n{}", dd_ep.as_matrix());

        // Compare the two tangent moduli (they should also equal the elastic stiffness)
        mat_approx_eq(dd_vm.matrix(), dd_ep.matrix(), 1e-15);
        mat_approx_eq(dd_vm.matrix(), ep.data_imp.as_ref().unwrap().args.dde.matrix(), 1e-15);

        // Set plane-strain strain increments such that the trial stress goes outside the yield surface
        let ee = YOUNG;
        let nu = POISSON;
        let nu2 = POISSON * POISSON;
        let z = 2.0 * Z_INI; // 2x the initial yield size so that the trial stress is outside the yield surface (it doesn't mean that we get twice the yield surface in the end because of plastic behavior)
        let dy = z * (1.0 - nu2) / (ee * f64::sqrt(1.0 - nu + nu2));
        let deps_x = dy * nu / (1.0 - nu);
        let deps_y = -dy;
        let mut delta_strain = Tensor2::new(mandel);
        delta_strain.vector_mut()[0] = deps_x;
        delta_strain.vector_mut()[1] = deps_y;

        // Array of states for plotting
        let mut states_vm = vec![state0.clone()];
        let mut states_ep = vec![state0.clone()];

        // Update the state using the von Mises model directly
        let mut state_vm = state0.clone();
        vm.update_stress(&mut state_vm, &delta_strain, 0, 0).unwrap();
        state_vm.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain); // eps += delta_eps
        states_vm.push(state_vm.clone());

        // Update the state using the Elastoplastic wrapper
        let mut state_ep = state0.clone();
        ep.update_stress(&mut state_ep, &delta_strain, 0, 0).unwrap();
        state_ep.strain.as_mut().unwrap().set_tensor(1.0, &delta_strain); // eps += delta_eps
        states_ep.push(state_ep.clone());

        // Compare the two updated states (they should be equal)
        vec_approx_eq(state_vm.stress.vector(), state_ep.stress.vector(), 1e-14);
        approx_eq(state_vm.int_vars[0], state_ep.int_vars[0], 1e-14);
        approx_eq(state_vm.lambda_alg, state_ep.lambda_alg, 1e-14);
        assert!(state_vm.lambda_alg > 0.0);

        // Calcualte the consistent tangent modulus at the updated state using the von Mises model
        let mut dd_vm = Tensor4::new(mandel);
        vm.stiffness(&mut dd_vm, &state_vm, 0, 0).unwrap();
        // println!("dd_vm =\n{}", dd_vm.as_matrix());

        // Calcualte the consistent tangent modulus at the updated state using the Elastoplastic wrapper
        let mut dd_ep = Tensor4::new(mandel);
        ep.stiffness(&mut dd_ep, &state_ep, 0, 0).unwrap();
        // println!("dd_ep =\n{}", dd_ep.as_matrix());

        // Compare the two tangent moduli
        mat_approx_eq(dd_vm.matrix(), dd_ep.matrix(), 1e-12);

        // plot
        if SAVE_FIGURE {
            let data_vm = PlotterData::from_states(&states_vm);
            let data_ep = PlotterData::from_states(&states_ep);
            let mut plotter = Plotter::new();
            plotter
                .add_2x2(&data_vm, false, |curve, _, _| {
                    curve.set_label("von Mises").set_marker_style("s").set_marker_void(true);
                })
                .unwrap();
            plotter
                .add_2x2(&data_ep, false, |curve, _, _| {
                    curve
                        .set_label("Elastoplastic")
                        .set_line_style(":")
                        .set_marker_style("o")
                        .set_marker_void(true);
                })
                .unwrap();
            let r0 = states_vm[0].int_vars[0] * SQRT_2_BY_3;
            let r1 = states_vm[1].int_vars[0] * SQRT_2_BY_3;
            plotter.set_oct_circle(r0, |_| {});
            plotter.set_oct_circle(r1, |canvas| {
                canvas.set_line_style("-");
            });
            plotter
                .save("/tmp/pmsim/material/test_implicit_consistent_modulus_matches_von_mises.svg")
                .unwrap();
        }
    }
}
