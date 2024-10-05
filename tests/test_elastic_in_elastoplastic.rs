use pmsim::base::{Idealization, StressStrain};
use pmsim::material::{Axis, Plotter, PlotterData, Settings, StressStrainTrait, VonMises};
use pmsim::material::{Elastoplastic, LinearElastic, LocalState};
use pmsim::util::elastic_increments_oct;
use pmsim::StrError;
use russell_lab::math::PI;
use russell_lab::vec_approx_eq;
use russell_tensor::{Mandel, Tensor2, SQRT_2_BY_3, SQRT_3_BY_2};

///////////////////////////////////////////////////////////////////////
//                                                                   //
// This tests verifies the elastic update calculated by the          //
// general Elastoplastic and von Mises models.                       //
//                                                                   //
///////////////////////////////////////////////////////////////////////

const FILE_STEM: &str = "test_elastic_in_elastoplastic";

const SAVE_FIGURE: bool = true;

const N_STEP: usize = 4;

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
    let mut box_elast: Box<dyn StressStrainTrait> = Box::new(elast);
    let mut box_direct: Box<dyn StressStrainTrait> = Box::new(direct);
    let mut box_general: Box<dyn StressStrainTrait> = Box::new(general);

    // constants
    let mandel = ideal.mandel();
    let n_int_val = box_direct.n_internal_values();
    assert_eq!(box_general.n_internal_values(), n_int_val);

    // initial states and increments
    let sig_m_0 = 1.0;
    let alpha_0 = PI / 3.0;
    let (stresses, strain_increments) = walk_on_oct_plane(young, poisson, z_ini, sig_m_0, alpha_0);

    // run test
    for i in 0..stresses.len() {
        // initial state
        let mut state_elast = LocalState::new(mandel, n_int_val);
        state_elast.stress.set_tensor(1.0, &stresses[i]);
        state_elast.enable_strain();
        let mut state_direct = state_elast.clone();
        let mut state_general = state_elast.clone();
        let mut state_general_full = state_elast.clone();

        // run test with sub-steps
        let depsilon = &strain_increments[i];
        let states_elast = update_with_steps(&mut state_elast, &mut box_elast, &depsilon, N_STEP)?;
        let states_direct = update_with_steps(&mut state_direct, &mut box_direct, &depsilon, N_STEP)?;
        let states_general = update_with_steps(&mut state_general, &mut box_general, &depsilon, N_STEP)?;

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
        for j in 0..N_STEP {
            let correct_stress = states_elast[j].stress.vector();
            vec_approx_eq(states_direct[j].stress.vector(), correct_stress, 1e-15);
            vec_approx_eq(states_general[j].stress.vector(), correct_stress, 1e-14);
        }

        // plot
        if SAVE_FIGURE {
            do_plot(i, &states_elast, &states_direct, &states_general, &general_full)?;
        }
    }
    Ok(())
}

// Updates stresses and strains using (sub)steps
fn update_with_steps(
    state: &mut LocalState,
    model: &mut Box<dyn StressStrainTrait>,
    depsilon_total: &Tensor2,
    n_step: usize,
) -> Result<Vec<LocalState>, StrError> {
    // check
    if n_step < 1 {
        return Err("n_step must be ≥ 1");
    }

    // initialize internal values
    model.initialize_internal_values(state)?;
    let mut states = vec![state.clone()];

    // update stress and strain
    let mut depsilon = Tensor2::new(state.stress.mandel());
    depsilon.update(1.0 / (n_step as f64), &depsilon_total);
    for _ in 0..n_step {
        model.update_stress(state, &depsilon)?;
        state.strain.as_mut().unwrap().update(1.0, &depsilon);
        states.push(state.clone())
    }
    Ok(states)
}

fn do_plot(
    index: usize,
    states_elast: &Vec<LocalState>,
    states_direct: &Vec<LocalState>,
    states_general: &Vec<LocalState>,
    general_full: &Elastoplastic,
) -> Result<(), StrError> {
    // constants
    let n = states_elast.len();
    let l = n - 1;
    let z_ini = states_general[0].internal_values[0];
    let z_fin = states_general[l].internal_values[0];

    // plotting data
    let mut data_elast = PlotterData::from_states(&states_elast);
    let mut data_direct = PlotterData::from_states(&states_direct);
    let mut data_general = PlotterData::from_states(&states_general);
    let tf = |i| (i as f64) / (l as f64);
    data_elast.set_time_and_yield(|i| Ok((tf(i), 0.0)))?;
    data_direct.set_time_and_yield(|i| Ok((tf(i), general_full.yield_function(&states_direct[i])?)))?;
    data_general.set_time_and_yield(|i| Ok((tf(i), general_full.yield_function(&states_general[i])?)))?;

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
            curve.set_label("direct").set_marker_style("o").set_marker_void(false);
        })
        .unwrap();
    plotter
        .add_3x2(&data_general, false, |curve, _, _| {
            curve.set_label("general").set_marker_style("s").set_marker_void(true);
        })
        .unwrap();
    let radius_0 = z_ini * SQRT_2_BY_3;
    let radius_1 = z_fin * SQRT_2_BY_3;
    plotter.set_oct_circle(radius_0, |_| {});
    plotter.set_oct_circle(radius_1, |canvas| {
        canvas.set_line_style("-");
    });
    plotter.set_extra(Axis::SigM(false), Axis::SigD(false), |plot| {
        plot.set_xrange(0.0, 2.0);
    });
    plotter.set_extra(Axis::EpsD(true), Axis::EpsV(true, false), |plot| {
        plot.set_yrange(-1.0, 1.0);
    });
    plotter.set_extra(Axis::SigM(false), Axis::EpsV(true, false), |plot| {
        plot.set_xrange(0.0, 2.0);
        plot.set_yrange(-1.0, 1.0);
    });
    plotter.set_figure_size(800.0, 1000.0);
    plotter.save(&format!("/tmp/pmsim/{}_{}.svg", FILE_STEM, index))
}

// Generates stresses and strain increments to "walk" on the octahedral plane along three directions
fn walk_on_oct_plane(young: f64, poisson: f64, z_ini: f64, sig_m_0: f64, alpha_0: f64) -> (Vec<Tensor2>, Vec<Tensor2>) {
    let r = z_ini * SQRT_2_BY_3;
    let rc = r * f64::cos(alpha_0);
    let rs = r * f64::sin(alpha_0);
    let mandel = Mandel::Symmetric2D;
    let mut initial_stresses = Vec::new();
    let mut strain_increments = Vec::new();
    for (oct_x_1, oct_y_1) in [
        (rc, -rs),  // bottom right
        (-rc, -rs), // bottom left
        (-rc, rs),  // top left
    ] {
        let alpha_1 = f64::atan2(oct_y_1, oct_x_1);
        let radius_1 = f64::sqrt(oct_x_1 * oct_x_1 + oct_y_1 * oct_y_1);
        let sig_d_1 = radius_1 * SQRT_3_BY_2;
        let (stress_0, _, _, d_strain) = elastic_increments_oct(
            young, poisson, sig_m_0, z_ini, alpha_0, sig_m_0, sig_d_1, alpha_1, mandel,
        );
        initial_stresses.push(stress_0);
        strain_increments.push(d_strain);
    }
    (initial_stresses, strain_increments)
}
