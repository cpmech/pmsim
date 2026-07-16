use crate::material::{LocalState, PlotterData, Settings, StressStrainTrait};
use super::data_exp::{Case, DataExp};
use crate::base::{Idealization, StressStrain};
use crate::StrError;
use gemlab::mesh::CellId;
use russell_tensor::{Tensor2, Tensor4};

/// Implements general elastoplasticity models using explicit stress update
pub struct ElastoplasticExp<'a> {
    /// Holds the data for the explicit stress update algorithm
    data: DataExp<'a>,

    /// Enables verbose mode
    pub verbose: bool,
}

impl<'a> ElastoplasticExp<'a> {
    /// Allocates a new instance
    pub fn new(ideal: &Idealization, param: &StressStrain, settings: &Settings) -> Result<Self, StrError> {
        Ok(ElastoplasticExp {
            data: DataExp::new(ideal, param, settings)?,
            verbose: false,
        })
    }

    /// Calculates the yield function f
    pub fn yield_function(&self, state: &LocalState) -> Result<f64, StrError> {
        self.data.args.model.calc_f(state)
    }

    /// Returns the stress-strain history during the intersection finding (e.g., for debugging)
    pub fn get_history_int(&self) -> Result<PlotterData, StrError> {
        match self.data.args.history_int.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns the stress-strain history during the elastic and elastoplastic update (e.g., for debugging)
    pub fn get_history_eep(&self) -> Result<PlotterData, StrError> {
        match self.data.args.history_eep.as_ref() {
            Some(h) => Ok(h.clone()),
            None => Err("history needs to be enabled (explicit update only)"),
        }
    }

    /// Returns the last Case analyzed by update_stress (for debugging)
    pub fn last_case(&self) -> Option<Case> {
        self.data.last_case
    }
}

impl<'a> StressStrainTrait for ElastoplasticExp<'a> {
    fn symmetric_stiffness(&self) -> bool {
        self.data.args.model.symmetric_stiffness()
    }

    fn n_int_vars(&self) -> usize {
        self.data.args.model.n_int_vars()
    }

    fn initialize_int_vars(&self, state: &mut LocalState) -> Result<(), StrError> {
        self.data.args.model.initialize_int_vars(state)
    }

    fn stiffness(
        &mut self,
        dd: &mut Tensor4,
        state: &LocalState,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        self.data.explicit_stiffness(dd, state)
    }

    fn update_stress(
        &mut self,
        state: &mut LocalState,
        delta_strain: &Tensor2,
        _cell_id: CellId,
        _gauss_id: usize,
    ) -> Result<(), StrError> {
        self.data.explicit_update_stress(state, delta_strain, self.verbose)
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Case, ElastoplasticExp};
    use crate::base::{Idealization, StressStrain};
    use crate::material::testing::{extract_von_mises_params, extract_von_mises_params_kg};
    use crate::material::{Axis, LocalState, Plotter, PlotterData, Settings, StressStrainTrait};
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
        let n_int_vars = model.n_int_vars();
        let mut state = LocalState::new(ideal.mandel(), n_int_vars);
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
        let param = StressStrain::sample_von_mises();
        let mut settings = Settings::new();
        settings.set_gp_explicit_update(true);
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 2.0);
        let (sig_m_1, sig_d_1) = (1.0, z_ini);
        let (sig_m_2, sig_d_2) = (2.0, 2.0 * z_ini);
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
                approx_eq(state.int_vars[0], z_ini, 1e-15);
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
                approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
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
        let (kk, gg, hh, z_ini) = extract_von_mises_params_kg(&param);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, 0.0, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, z_ini + 9.0, PI / 3.0);
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
        let deps_d_e = z_ini / (3.0 * gg);
        let deps_d_ep = deps_d - deps_d_e;
        let correct_sig_m = kk * deps_v;
        let correct_sig_d = z_ini + 3.0 * gg * hh * deps_d_ep / (3.0 * gg + hh);
        approx_eq(sig_m, correct_sig_m, 1e-14);
        approx_eq(sig_d, correct_sig_d, 1e-13);
        approx_eq(state.int_vars[0], correct_sig_d, 1e-13);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (0.0, 0.99999999999999999);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0);
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
        approx_eq(state.int_vars[0], z_ini, 1e-15);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1.0, 0.8);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -2.0 * PI / 3.0);
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
        approx_eq(state.int_vars[0], z_ini, 1e-15);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
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
        approx_eq(state.int_vars[0], z_ini, 1e-15);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let drift = 0.0;
        let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
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
        approx_eq(state.int_vars[0], z_ini, 1e-15);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1.0, 2.5);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, -PI / 3.0);
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
        let (_, _, _, z_ini) = extract_von_mises_params_kg(&param);
        let (drift, mz) = (1e-14, 2.0);
        let (sig_m_0, sig_d_0, alpha_0) = (0.0, z_ini + drift, PI / 3.0);
        let (sig_m_1, sig_d_1, alpha_1) = (2.0, mz * z_ini, 0.0);
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
