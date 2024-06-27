use super::{StressStrainModelName, StressStrainPlot, StressStrainState};
use plotpy::{Canvas, Curve, Legend, RayEndpoint};
use russell_lab::math::SQRT_2_BY_3;

pub struct GraphElastoplastic {
    pub n_markers: usize,
    pub out_dir: String,
    pub color_yield_f_ini: String,
    pub color_yield_f_fin: String,
}

impl GraphElastoplastic {
    pub fn new() -> GraphElastoplastic {
        GraphElastoplastic {
            n_markers: 4,
            out_dir: "/tmp/pmsim/material".to_string(),
            color_yield_f_ini: "#CD54DF".to_string(),
            color_yield_f_fin: "#B12AC4".to_string(),
        }
    }

    pub fn standard_vs_general(
        &self,
        filename_key: &str,
        model_name: StressStrainModelName,
        std_states: &Vec<StressStrainState>,
        gen_states: &Vec<StressStrainState>,
        gen_history_e: Option<&Vec<StressStrainState>>,
        gen_history_ep: Option<&Vec<StressStrainState>>,
    ) {
        let mut ssp = StressStrainPlot::new();
        if let Some(history_e) = gen_history_e {
            let mark_every = history_e.len() / self.n_markers;
            ssp.draw_2x2_mosaic_struct(history_e, |curve, _, _| {
                curve
                    .set_marker_style("s")
                    .set_line_style(":")
                    .set_label("history(e)")
                    .set_marker_every(mark_every);
            });
        }
        if let Some(history_ep) = gen_history_ep {
            let mark_every = history_ep.len() / self.n_markers;
            ssp.draw_2x2_mosaic_struct(history_ep, |curve, _, _| {
                curve
                    .set_marker_style("^")
                    .set_line_style(":")
                    .set_label("history(ep)")
                    .set_marker_every(mark_every);
            });
        }
        ssp.draw_2x2_mosaic_struct(gen_states, |curve, _row, _col| {
            curve
                .set_marker_style("o")
                .set_marker_void(true)
                .set_line_style(":")
                .set_marker_size(10.0)
                .set_label("general");
        });
        ssp.draw_2x2_mosaic_struct(std_states, |curve, _row, _col| {
            curve.set_marker_style("o").set_label("standard");
        });
        let mut legend = Legend::new();
        legend.set_outside(true).set_num_col(2);
        let filepath = &format!("{}/{}.svg", self.out_dir, filename_key);
        ssp.save_2x2_mosaic_struct(filepath, |plot, row, col, before| {
            if before {
                if col == 0 {
                    match model_name {
                        StressStrainModelName::VonMises => {
                            let z_ini = std_states.first().unwrap().internal_values[0];
                            let z_fin = std_states.last().unwrap().internal_values[0];
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
                            let z_ini = std_states.first().unwrap().internal_values[0];
                            let z_fin = std_states.last().unwrap().internal_values[0];
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
                    plot.set_cross(0.0, 0.0, "grey", "-", 2.0);
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
