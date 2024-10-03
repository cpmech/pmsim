use plotpy::Text;
use pmsim::base::{Idealization, StressStrain};
use pmsim::material::{Axis, Plotter, PlotterData, Settings, StressStrainTrait, VonMises};
use pmsim::material::{Elastoplastic, LinearElastic, LocalState};
use pmsim::StrError;
use russell_lab::math::PI;
use russell_tensor::{Tensor2, SQRT_2_BY_3, SQRT_3};

const FILE_STEM: &str = "test_elastic_in_elastoplastic";

fn run(
    state: &mut LocalState,
    mut model: Box<dyn StressStrainTrait>,
    depsilon_total: &Tensor2,
    n_step: usize,
) -> Result<PlotterData, StrError> {
    // check
    if n_step < 1 {
        return Err("n_step must be ≥ 1");
    }

    // initialize internal values
    model.initialize_internal_values(state)?;

    // initialize plotting data
    let mut data = PlotterData::new();
    data.push(
        &state.stress,
        Some(state.strain.as_ref().unwrap()),
        Some(0.0),
        Some(0.0),
    );

    // sub-steps
    let mut depsilon = Tensor2::new(state.stress.mandel());
    depsilon.update(1.0 / (n_step as f64), &depsilon_total);
    for _ in 0..n_step {
        // update stress and strain
        model.update_stress(state, &depsilon)?;
        state.strain.as_mut().unwrap().update(1.0, &depsilon);

        // update plotting data
        data.push(
            &state.stress,
            Some(state.strain.as_ref().unwrap()),
            Some(0.0),
            Some(1.0),
        );
    }
    Ok(data)
}

#[test]
fn test_elastic_in_elastoplastic() -> Result<(), StrError> {
    // parameters
    let young = 1500.0;
    let poisson = 0.25;
    let hh = 800.0;
    let z_ini = 9.0;
    let param_el = StressStrain::LinearElastic { young, poisson };
    let param_vm = StressStrain::VonMises {
        young,
        poisson,
        hh,
        z_ini,
    };

    // models
    let ndim = 2;
    let ideal = Idealization::new(ndim);
    let mut settings = Settings::new();
    settings.set_gp_save_history(true);
    let elast = LinearElastic::new(&ideal, &param_el, &settings)?;
    let direct = VonMises::new(&ideal, &param_vm, &settings)?;
    let general = Elastoplastic::new(&ideal, &param_vm, &settings)?;
    let mut general_full = Elastoplastic::new(&ideal, &param_vm, &settings)?;
    let box_elast: Box<dyn StressStrainTrait> = Box::new(elast);
    let box_direct: Box<dyn StressStrainTrait> = Box::new(direct);
    let box_general: Box<dyn StressStrainTrait> = Box::new(general);

    // constants
    let mandel = ideal.mandel();
    let two_dim = ideal.two_dim;
    let n_int_val = box_direct.n_internal_values();
    assert_eq!(box_general.n_internal_values(), n_int_val);

    // initial state
    let sig_m_0 = 1.0;
    let sig_d_0 = z_ini;
    let alpha_0 = PI / 3.0;
    let dist_0 = sig_m_0 * SQRT_3;
    let radius_0 = sig_d_0 * SQRT_2_BY_3;
    let mut state_elast = LocalState::new(mandel, n_int_val);
    state_elast.stress = Tensor2::new_from_octahedral_alpha(dist_0, radius_0, alpha_0, two_dim)?;
    state_elast.enable_strain();
    let mut state_direct = state_elast.clone();
    let mut state_general = state_elast.clone();
    let mut state_general_full = state_elast.clone();

    // total strain increment
    let depsilon = Tensor2::from_matrix(
        &[
            [-0.008660254037844387, 0.0, 0.0],
            [0.0, 0.004330127018922193, 0.0],
            [0.0, 0.0, 0.004330127018922193],
        ],
        mandel,
    )?;

    // run test with sub-steps
    let n_step = 7;
    let data_elast = run(&mut state_elast, box_elast, &depsilon, n_step)?;
    let data_direct = run(&mut state_direct, box_direct, &depsilon, n_step)?;
    let data_general = run(&mut state_general, box_general, &depsilon, n_step)?;

    // run general with full step
    let mut data_general_full = PlotterData::new();
    data_general_full.push(
        &state_general_full.stress,
        Some(state_general_full.strain.as_ref().unwrap()),
        Some(general_full.yield_function(&state_general_full)?),
        Some(0.0),
    );
    general_full.initialize_internal_values(&mut state_general_full)?;
    general_full.update_stress(&mut state_general_full, &depsilon)?;
    data_general_full.push(
        &state_general_full.stress,
        Some(state_general_full.strain.as_ref().unwrap()),
        Some(general_full.yield_function(&state_general_full)?),
        Some(1.0),
    );

    // check
    // TODO

    // plot results
    let mut plotter = Plotter::new();
    plotter
        .set_tab_leg_ncol(2)
        .set_layout_selected_3x2(Axis::Time, Axis::Yield);
    plotter
        .add_3x2(&data_elast, false, |curve, _, _| {
            curve
                .set_label("elastic")
                .set_line_color("#01710e")
                .set_marker_size(12.0)
                .set_marker_void(true)
                .set_marker_style("h");
        })
        .unwrap();
    let history_int = general_full.get_history_int().unwrap();
    let history_eep = general_full.get_history_eep().unwrap();
    plotter
        .add_3x2(&history_int, false, |curve, _, _| {
            curve
                .set_label("history(int)")
                .set_line_color("purple")
                .set_line_style("-");
        })
        .unwrap();
    plotter
        .add_3x2(&history_eep, false, |curve, _, _| {
            curve
                .set_label("history(e-ep)")
                .set_line_color("#7a7a7a")
                .set_line_style("--")
                .set_marker_style(".")
                .set_marker_every(2);
        })
        .unwrap();
    plotter
        .add_3x2(&data_direct, false, |curve, _, _| {
            curve.set_label("direct").set_marker_style("o").set_marker_void(true);
        })
        .unwrap();
    plotter
        .add_3x2(&data_general, false, |curve, _, _| {
            curve.set_label("general").set_marker_style("s").set_marker_void(true);
        })
        .unwrap();
    let radius_0 = z_ini * SQRT_2_BY_3;
    let radius_1 = state_general.internal_values[0] * SQRT_2_BY_3;
    plotter.set_oct_circle(radius_0, |_| {});
    plotter.set_oct_circle(radius_1, |canvas| {
        canvas.set_line_style("-");
    });
    plotter.set_extra(Axis::SigM(false), Axis::SigD(false), |plot| {
        plot.set_xrange(0.0, 2.0);
    });
    plotter.set_extra(Axis::SigM(false), Axis::EpsV(true, false), |plot| {
        plot.set_xrange(0.0, 2.0);
    });
    plotter.set_extra(Axis::OctX, Axis::OctY, move |_plot| {
        // let mut text = get_text_label();
        // plot.add(&text);
    });
    plotter.set_extra(Axis::Time, Axis::Yield, move |_plot| {
        // let mut text = get_text_label();
        // plot.add(&text);
    });
    plotter.set_figure_size(800.0, 1000.0);
    plotter.save(&format!("/tmp/pmsim/{}.svg", FILE_STEM)).unwrap();

    Ok(())
}

// Returns Text for labels in plots
fn _get_text_label() -> Text {
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

/*
let (sig_m_0, sig_d_0, alpha_0) = (1.0, z_ini + drift, PI / 3.0);
let sig_m_1 = sig_m_0;
let radius_0 = sig_d_0 * SQRT_2_BY_3;
let (oct_x_1, oct_y_1) = (radius_0 * f64::cos(alpha_0), -radius_0 * f64::sin(alpha_0));
let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
let sig_d_1 = radius_1 * SQRT_3_BY_2;
*/
