use super::{callback_history_e, callback_history_ep, callback_intersect, callback_ode_e, callback_ode_ep};
use super::{Args, Case};
use super::{CHEBYSHEV_TOL, HISTORY_N_OUT, PSEUDO_TIME_TOL};
use crate::base::{Idealization, StressStrain};
use crate::material::{LocalState, PlotterData, Settings, TraitStressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_lab::{InterpChebyshev, RootFinder, Vector};
use russell_ode::{OdeSolver, Output, Params, System};
use russell_tensor::{t2_ddot_t4_ddot_t2, Tensor2, Tensor4};

/// Implements general elastoplasticity models using explicit stress update
pub struct ElastoplasticExp<'a> {
    /// Holds the arguments for the explicit stress update algorithm
    args: Args,

    /// Holds the solver for finding the yield surface intersection
    ode_intersection: OdeSolver<'a, Args>,

    /// Holds the solver for the elastic update
    ode_elastic: OdeSolver<'a, Args>,

    /// Holds the solver for the elastoplastic update
    ode_elastoplastic: OdeSolver<'a, Args>,

    /// Holds the ODE vector of unknowns for elastic case
    ode_y_e: Vector,

    /// Holds the ODE vector of unknowns for elastoplastic case
    ode_y_ep: Vector,

    /// Holds the output during the intersection finding
    out_intersection: Output<'a, Args>,

    /// Holds the output during the elastic path
    out_history_el: Output<'a, Args>,

    /// Holds the output during the elastoplastic path
    out_history_ep: Output<'a, Args>,

    /// Holds the interpolant for finding the yield surface intersection
    interpolant: InterpChebyshev,

    /// Solver for the intersection finding algorithm
    root_finder: RootFinder,

    /// Enables recording stress-strain history
    save_history: bool,

    /// Holds the last Case analyzed by update_stress (for debugging)
    last_case: Option<Case>,

    /// Enables verbose mode
    verbose: bool,
}

impl<'a> ElastoplasticExp<'a> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        // Allocate the interpolant
        let interp_nn_max = settings.gp_interp_nn_max();
        let interpolant = InterpChebyshev::new(interp_nn_max, 0.0, 1.0).unwrap();

        // Allocate the Chebyshev points and the interior t_out values
        let chebyshev_points = InterpChebyshev::points(interp_nn_max);
        let interp_npoint = chebyshev_points.dim();
        let mut interior_t_out = vec![0.0; interp_npoint - 2];
        let xx_interior = &chebyshev_points.as_data()[1..(interp_npoint - 1)];
        xx_interior.into_iter().enumerate().for_each(|(i, x)| {
            interior_t_out[i] = (1.0 + x) / 2.0;
        });

        // Allocate the arguments for the explicit stress update algorithm
        let args = Args::new(ideal, param, settings, interp_npoint)?;

        // Allocate the ODE systems
        let ode_system_e = System::new(args.ndim_e, callback_ode_e);
        let ode_system_ep = System::new(args.ndim_ep, callback_ode_ep);

        // Allocate the ODE solvers
        let ode_param = Params::new(settings.gp_ode_method());
        let ode_intersection = OdeSolver::new(ode_param, ode_system_e.clone()).unwrap();
        let ode_elastic = OdeSolver::new(ode_param, ode_system_e).unwrap();
        let ode_elastoplastic = OdeSolver::new(ode_param, ode_system_ep).unwrap();

        // Allocate the output for the intersection finding
        let mut out_intersection = Output::new();
        out_intersection
            .set_dense_x_out(&interior_t_out)
            .unwrap()
            .set_dense_callback(callback_intersect);

        // Allocate the output for the elastic and elastoplastic paths
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

        // Allocate the ODE vectors of unknowns
        let ode_y_e = Vector::new(args.ndim_e);
        let ode_y_ep = Vector::new(args.ndim_ep);
        let root_finder = RootFinder::new();

        // Allocate the instance
        Ok(ElastoplasticExp {
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
            verbose: settings.gp_verbose(),
        })
    }

    /// Calculates the yield function f
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        self.args.model.calc_f(state)
    }

    /// Returns the stress-strain history during the intersection finding (e.g., for debugging)
    pub fn get_history_int(&self) -> Result<PlotterData, StrError> {
        match self.args.history_int.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub fn get_history_eep(&self) -> Result<PlotterData, StrError> {
        match self.args.history_eep.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns the last Case analyzed by update_stress (for debugging)
    pub fn last_case(&self) -> Option<Case> {
        self.last_case
    }

    /// Returns true if the trial stress path leads to the inside of the yield surface
    fn going_inside(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<bool, StrError> {
        self.args.model.calc_fs(&mut self.args.fs, state)?;
        self.args.model.calc_dde(&mut self.args.dde, state)?;
        let indicator = t2_ddot_t4_ddot_t2(&self.args.fs, &self.args.dde, delta_strain);
        Ok(indicator < 0.0)
    }

    /// Performs the intersection finding algorithm
    fn intersection_finding(&mut self, state: &LocalState, inside: bool) -> Result<(Option<f64>, f64), StrError> {
        self.args.state.z_set.set_vector(state.z_set.as_data());
        self.ode_y_e.set_vector(state.stress.vector().as_data());
        self.ode_intersection.solve(
            &mut self.ode_y_e,
            0.0,
            1.0,
            None,
            &mut self.args,
            Some(&mut self.out_intersection),
        )?;
        assert_eq!(self.args.yf_count, self.args.yf_values.dim());
        self.interpolant
            .adapt_data(CHEBYSHEV_TOL, self.args.yf_values.as_data())?;
        let roots = self.root_finder.chebyshev(&self.interpolant)?;
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
        let yf_trial = self.args.yf_values[self.args.yf_count - 1];
        Ok((t_int, yf_trial))
    }

    /// Selects the yield surface crossing case
    fn select_case(&mut self, state: &LocalState, delta_strain: &Tensor2) -> Result<Case, StrError> {
        let yf_initial = self.args.model.calc_f(state)?;
        if yf_initial < 0.0 {
            let (t_intersection, yf_trial) = self.intersection_finding(state, true)?;
            match t_intersection {
                Some(t_int) => {
                    if t_int <= PSEUDO_TIME_TOL {
                        Ok(Case::BP)
                    } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                        Ok(Case::AE)
                    } else {
                        Ok(Case::AXB(t_int))
                    }
                }
                None => {
                    assert!(yf_trial <= 0.0);
                    Ok(Case::AE)
                }
            }
        } else {
            if self.going_inside(state, delta_strain)? {
                let (t_intersection, yf_trial) = self.intersection_finding(state, false)?;
                match t_intersection {
                    Some(t_int) => {
                        if t_int <= PSEUDO_TIME_TOL {
                            Ok(Case::BP)
                        } else if t_int >= 1.0 - PSEUDO_TIME_TOL {
                            Ok(Case::BE)
                        } else {
                            Ok(Case::BXP(t_int))
                        }
                    }
                    None => {
                        assert!(yf_trial <= 0.0);
                        Ok(Case::BE)
                    }
                }
            } else {
                Ok(Case::BP)
            }
        }
    }
}

impl<'a> TraitStressStrain for ElastoplasticExp<'a> {
    /// Returns whether this model has symmetric stiffness matrix or not
    fn symmetric_stiffness(&self) -> bool {
        self.args.model.symmetric_stiffness()
    }

    /// Returns the number of internal variables
    fn nz(&self) -> usize {
        self.args.model.nz()
    }

    /// Initializes the internal variables for the initial stress state
    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.args.model.initialize_int_vars(state)
    }

    /// Returns an error because the stiffness is not available for explicit update
    fn stiffness(
        &mut self,
        _dd: &mut Tensor4,
        _state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        Err("stiffness is not available for explicit update")
    }

    /// Updates the stress tensor given the strain increment tensor using the explicit method
    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        {
            self.args.del_eps.set_tensor(1.0, delta_strain);
            if self.save_history {
                match state.strain.as_ref() {
                    Some(strain) => self.args.state.strain = Some(strain.clone()),
                    None => {
                        return Err("state must have strain enabled");
                    }
                }
                self.args.history_int = Some(PlotterData::new());
                self.args.history_eep = Some(PlotterData::new());
            }
        }

        let case = self.select_case(state, delta_strain)?;

        match case {
            Case::AE | Case::BE => {
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                state.elastic = true;
            }
            Case::AXB(t_int) | Case::BXP(t_int) => {
                self.ode_y_e.set_vector(state.stress.vector().as_data());
                self.ode_elastic.solve(
                    &mut self.ode_y_e,
                    0.0,
                    t_int,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_el),
                )?;
                state.stress.vector_mut().set_vector(self.ode_y_e.as_data());
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.z_set.as_data());
                self.ode_elastoplastic.solve(
                    &mut self.ode_y_ep,
                    t_int,
                    1.0,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_ep),
                )?;
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.z_set.as_mut_data());
                state.elastic = false;
            }
            Case::BP => {
                self.ode_y_ep
                    .join2(state.stress.vector().as_data(), state.z_set.as_data());
                self.ode_elastoplastic.solve(
                    &mut self.ode_y_ep,
                    0.0,
                    1.0,
                    None,
                    &mut self.args,
                    Some(&mut self.out_history_ep),
                )?;
                self.ode_y_ep
                    .split2(state.stress.vector_mut().as_mut_data(), state.z_set.as_mut_data());
                state.elastic = false;
            }
        }

        if self.verbose {
            println!("👉 {:?}", case);
        }

        self.last_case = Some(case);
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Case, ElastoplasticExp};
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::{extract_von_mises_params, extract_von_mises_params_kg};
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, TraitStressStrain};
    use plotpy::Text;
    use russell_lab::{approx_eq, math::PI};
    use russell_tensor::{t2_add, t4_ddot_t2, LinElasticity, Tensor2, Tensor4, SQRT_2_BY_3, SQRT_3, SQRT_3_BY_2};
    use std::collections::HashMap;

    const VERBOSE: bool = true;
    const SAVE_FIGURE: bool = false;

    fn case_to_keys(case: Case) -> Vec<&'static str> {
        match case {
            Case::AE => vec!["A", "E"],
            Case::AXB(..) => vec!["A", "X", "B"],
            Case::BE => vec!["B", "E"],
            Case::BXP(..) => vec!["B", "X", "P"],
            Case::BP => vec!["B", "P"],
        }
    }

    fn gen_ini_state_von_mises(
        ideal: &Idealization,
        model: &ElastoplasticExp,
        p: f64,
        q: f64,
        alpha: f64,
    ) -> LocalState {
        let distance = p * SQRT_3;
        let radius = q * SQRT_2_BY_3;
        let nz = model.nz();
        let mut state = LocalState::new(ideal.mandel(), nz);
        state.stress = Tensor2::new_from_octahedral_alpha(distance, radius, alpha, ideal.two_dim).unwrap();
        model.initialize_int_vars(&mut state).unwrap();
        state.enable_strain();
        state
    }

    fn update_with_von_mises(
        param: &StressStrain,
        model: &mut ElastoplasticExp,
        state: &mut LocalState,
        p_el: f64,
        q_el: f64,
        alpha_el: f64,
    ) -> (f64, f64) {
        let distance = p_el * SQRT_3;
        let radius = q_el * SQRT_2_BY_3;
        let mandel = state.stress.mandel();
        let two_dim = mandel.two_dim();
        let stress_fin = Tensor2::new_from_octahedral_alpha(distance, radius, alpha_el, two_dim).unwrap();
        let mut dsigma = Tensor2::new(mandel);
        t2_add(&mut dsigma, 1.0, &stress_fin, -1.0, &state.stress);
        let (young, poisson, _, _) = extract_von_mises_params(param);
        let elast = LinElasticity::new(young, poisson, two_dim, false);
        let mut cc = Tensor4::new(mandel);
        elast.calc_compliance(&mut cc).unwrap();
        let mut depsilon = Tensor2::new(mandel);
        t4_ddot_t2(&mut depsilon, 1.0, &cc, &dsigma);
        model.update_stress(state, &depsilon, 0, 0).unwrap();
        state.strain.as_mut().unwrap().update(1.0, &depsilon);
        (depsilon.invariant_eps_v(), depsilon.invariant_eps_d())
    }

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
                let f = s.stress.invariant_q() - s.z_set[0];
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
                let radius_0 = states[0].z_set[0] * SQRT_2_BY_3;
                let radius_1 = states[p].z_set[0] * SQRT_2_BY_3;
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

    fn do_plot_b(
        file_stem: &str,
        model: &ElastoplasticExp,
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
        let radius_0 = states[0].z_set[0] * SQRT_2_BY_3;
        let radius_1 = states[p].z_set[0] * SQRT_2_BY_3;
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
        let param = StressStrain::sample_von_mises();
        let mut settings = Settings::new();
        settings.set_gp_explicit_update(true);
        let (kk, gg, hh, kappa_ini) = extract_von_mises_params_kg(&param);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 2.0);
        let (sig_m_1, sig_d_1) = (1.0, kappa_ini);
        let (sig_m_2, sig_d_2) = (2.0, 2.0 * kappa_ini);
        let mut data_2d = HashMap::new();
        for ndim in [2, 3] {
            for lode_int in [-1, 0, 1] {
                let lode = lode_int as f64;
                let alpha = PI / 2.0 - f64::acos(lode) / 3.0;
                let alpha_deg = alpha * 180.0 / PI;
                if VERBOSE {
                    println!("\nndim = {}, lode = {}, alpha = {}°", ndim, lode, alpha_deg);
                }
                let ideal = Idealization::new(ndim);
                let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
                model.verbose = VERBOSE;
                let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
                if ndim == 2 {
                    data_2d.insert(lode_int, vec![state.clone()]);
                }
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha);
                let sig_m_1 = state.stress.invariant_p();
                let sig_d_1 = state.stress.invariant_q();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }
                let correct_sig_m = sig_m_0 + kk * deps_v;
                let correct_sig_d = sig_d_0 + 3.0 * gg * deps_d;
                approx_eq(sig_m_1, correct_sig_m, 1e-14);
                approx_eq(sig_d_1, correct_sig_d, 1e-13);
                approx_eq(state.z_set[0], kappa_ini, 1e-15);
                assert_eq!(state.elastic, true);
                let case = model.last_case().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, ["A", "E"]);
                let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_2, sig_d_2, alpha);
                let sig_m_2 = state.stress.invariant_p();
                let sig_d_2 = state.stress.invariant_q();
                if ndim == 2 {
                    data_2d.get_mut(&lode_int).unwrap().push(state.clone());
                }
                let correct_sig_m = sig_m_1 + kk * deps_v;
                let correct_sig_d = sig_d_1 + 3.0 * gg * hh * deps_d / (3.0 * gg + hh);
                approx_eq(sig_m_2, correct_sig_m, 1e-14);
                approx_eq(sig_d_2, correct_sig_d, 1e-13);
                approx_eq(state.z_set[0], correct_sig_d, 1e-13);
                assert_eq!(state.elastic, false);
                let case = model.last_case().unwrap();
                let keys = case_to_keys(case);
                assert_eq!(keys, &["B", "P"]);
            }
        }
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
        let param = StressStrain::sample_von_mises();
        let (kk, gg, hh, kappa_ini) = extract_von_mises_params_kg(&param);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, kappa_ini + 9.0, PI / 3.0);
        let mut settings = Settings::new();
        settings.set_gp_explicit_update(true).set_gp_save_history(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        let (deps_v, deps_d) = update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());
        let deps_d_e = kappa_ini / (3.0 * gg);
        let deps_d_ep = deps_d - deps_d_e;
        let correct_sig_m = kk * deps_v;
        let correct_sig_d = kappa_ini + 3.0 * gg * hh * deps_d_ep / (3.0 * gg + hh);
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-13);
        approx_eq(state.z_set[0], correct_sig_d, 1e-13);
        assert_eq!(state.elastic, false);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["A", "X", "B"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (0.0, 0.99999999999999999);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, kappa_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * kappa_ini, -2.0 * PI / 3.0);
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.z_set[0], kappa_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "E"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1.0, 0.8);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, kappa_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * kappa_ini, -2.0 * PI / 3.0);
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.z_set[0], kappa_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "E"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, kappa_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (radius_0 * f64::cos(alpha_0), -radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.z_set[0], kappa_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "E"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, kappa_ini + drift, PI / 3.0);
        let sig_m_1 = sig_m_0;
        let radius_0 = sig_d_0 * SQRT_2_BY_3;
        let (oct_x_1, oct_y_1) = (-radius_0 * f64::cos(alpha_0), radius_0 * f64::sin(alpha_0));
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        let sig_m = state.stress.invariant_p();
        let sig_d = state.stress.invariant_q();
        states.push(state.clone());
        approx_eq(sig_m, sig_m_1, 1e-14);
        approx_eq(sig_d, sig_d_1, 1e-13);
        approx_eq(state.z_set[0], kappa_ini, 1e-15);
        assert_eq!(state.elastic, true);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "E"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1.0, 2.5);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, kappa_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * kappa_ini, -PI / 3.0);
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());
        assert_eq!(state.elastic, false);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "X", "P"]);
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
        let param = StressStrain::sample_von_mises();
        let (_, _, _, kappa_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1e-14, 2.0);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, kappa_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * kappa_ini, 0.0);
        let mut settings = Settings::new();
        settings
            .set_gp_explicit_update(true)
            .set_gp_save_history(true)
            .set_gp_allow_initial_drift(true);
        let ndim = 2;
        let ideal = Idealization::new(ndim);
        let mut model = ElastoplasticExp::new(&ideal, &param, &settings).unwrap();
        model.verbose = VERBOSE;
        let mut state = gen_ini_state_von_mises(&ideal, &model, sig_m_0, sig_d_0, alpha_0);
        let mut states = vec![state.clone()];
        update_with_von_mises(&param, &mut model, &mut state, sig_m_1, sig_d_1, alpha_1);
        states.push(state.clone());
        assert_eq!(state.elastic, false);
        let keys = case_to_keys(model.last_case().unwrap());
        assert_eq!(keys, &["B", "P"]);
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
}
