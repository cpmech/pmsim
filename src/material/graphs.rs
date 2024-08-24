use super::{LocalStateOld, StressStrainModelName, StressStrainPlot};
use plotpy::{Canvas, Curve, Legend, RayEndpoint};
use russell_lab::math::SQRT_2_BY_3;

pub struct GraphElastoplastic {
    pub n_markers: usize,
    pub out_dir: String,
    pub color_yield_f_ini: String,
    pub color_yield_f_fin: String,
    pub color_std: String,
    pub color_gen: String,
    pub color_gen_history_e: String,
    pub color_gen_history_ep: String,
    pub marker_style_history_e: String,
    pub marker_style_history_ep: String,
    pub marker_size_gen: f64,
    pub std_label: String,
    pub gen_label: String,
}

impl GraphElastoplastic {
    pub fn new() -> GraphElastoplastic {
        GraphElastoplastic {
            n_markers: 4,
            out_dir: "/tmp/pmsim/material".to_string(),
            color_yield_f_ini: "#FFA000".to_string(),
            color_yield_f_fin: "#EF6C00".to_string(),
            color_std: "#138D75".to_string(),
            color_gen: "#C62828".to_string(),
            color_gen_history_e: "#1E88E5".to_string(),
            color_gen_history_ep: "#1E88E5".to_string(),
            marker_style_history_e: "x".to_string(),
            marker_style_history_ep: "*".to_string(),
            marker_size_gen: 12.0,
            std_label: "return-map".to_string(),
            gen_label: "rk-general".to_string(),
        }
    }

    pub fn standard_vs_general(
        &self,
        filename_key: &str,
        model_name: StressStrainModelName,
        std_states: &Vec<LocalStateOld>,
        gen_states: &Vec<LocalStateOld>,
        gen_history_e: Option<&Vec<Vec<LocalStateOld>>>,
        gen_history_ep: Option<&Vec<Vec<LocalStateOld>>>,
    ) {
        let mut ssp = StressStrainPlot::new();
        ssp.no_grid = true;
        // general: history elastic
        if let Some(histories_e) = gen_history_e {
            let mut with_label = true;
            for history_e in histories_e {
                if history_e.len() > 0 {
                    let mark_every = history_e.len() / self.n_markers;
                    ssp.draw_2x2_mosaic_struct(history_e, |curve, _, _| {
                        curve
                            .set_line_color(&self.color_gen_history_e)
                            .set_marker_style(&self.marker_style_history_e)
                            .set_line_style("-")
                            .set_marker_every(mark_every);
                        if with_label {
                            curve.set_label("history(e)");
                        }
                    });
                    with_label = false;
                }
            }
        }
        // general: history elastoplastic
        if let Some(histories_ep) = gen_history_ep {
            let mut with_label = true;
            for history_ep in histories_ep {
                if history_ep.len() > 0 {
                    let mark_every = history_ep.len() / self.n_markers;
                    ssp.draw_2x2_mosaic_struct(history_ep, |curve, _, _| {
                        curve
                            .set_line_color(&self.color_gen_history_ep)
                            .set_marker_style(&self.marker_style_history_ep)
                            .set_line_style("-")
                            .set_marker_every(mark_every);
                        if with_label {
                            curve.set_label("history(ep)");
                        }
                    });
                    with_label = false;
                }
            }
        }
        // standard
        ssp.draw_2x2_mosaic_struct(std_states, |curve, _row, _col| {
            curve
                .set_line_color(&self.color_std)
                .set_line_style(":")
                .set_marker_style("o")
                .set_label(&self.std_label);
        });
        // general
        ssp.draw_2x2_mosaic_struct(gen_states, |curve, _row, _col| {
            curve
                .set_line_color(&self.color_gen)
                .set_line_style("None")
                .set_marker_style("o")
                .set_marker_void(true)
                .set_marker_size(self.marker_size_gen)
                .set_label(&self.gen_label);
        });
        // save figure
        let mut legend = Legend::new();
        legend.set_outside(true).set_num_col(2);
        let filepath = &format!("{}/{}.svg", self.out_dir, filename_key);
        ssp.save_2x2_mosaic_struct(filepath, |plot, row, col, before| {
            if before {
                if col == 0 {
                    match model_name {
                        StressStrainModelName::VonMises => {
                            let z_ini = gen_states.first().unwrap().internal_values[0];
                            let z_fin = gen_states.last().unwrap().internal_values[0];
                            let mut surf = Curve::new();
                            surf.set_line_color(&self.color_yield_f_ini)
                                .set_line_style("--")
                                .draw_ray(0.0, z_ini, RayEndpoint::Horizontal);
                            surf.set_line_color(&self.color_yield_f_fin)
                                .set_line_style("-")
                                .set_line_width(1.5)
                                .draw_ray(0.0, z_fin, RayEndpoint::Horizontal);
                            plot.add(&surf);
                        }
                        _ => (),
                    }
                }
                if row == 0 && col == 1 {
                    match model_name {
                        StressStrainModelName::VonMises => {
                            let z_ini = gen_states.first().unwrap().internal_values[0];
                            let z_fin = gen_states.last().unwrap().internal_values[0];
                            let mut surf = Canvas::new();
                            surf.set_face_color("None");
                            surf.set_edge_color(&self.color_yield_f_ini)
                                .set_line_style("--")
                                .draw_circle(0.0, 0.0, z_ini * SQRT_2_BY_3);
                            surf.set_edge_color(&self.color_yield_f_fin)
                                .set_line_style("-")
                                .draw_circle(0.0, 0.0, z_fin * SQRT_2_BY_3);
                            plot.add(&surf);
                        }
                        _ => (),
                    }
                }
                if row == 1 && col == 1 {
                    plot.set_cross(0.0, 0.0, "grey", "-", 1.3);
                }
            } else {
                if row == 1 && col == 1 {
                    legend.draw();
                    plot.add(&legend);
                }
            }
        })
        .unwrap();
    }
}
