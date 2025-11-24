use russell_lab::math::{elliptic_f, PI};

const COMPLETE: f64 = PI / 2.0;

/// Confined seepage flow underneath a dam using the Polubarinova-Kochina Solution
pub struct ConfinedFlow {}

impl ConfinedFlow {
    pub fn polubarinova_kochina_solution(s_by_tt: f64, b_by_tt: f64) -> f64 {
        if s_by_tt < 0.0 || s_by_tt > 1.0 {
            panic!("s/T must be between 0 and 1. {} is invalid", s_by_tt);
        }
        if b_by_tt < 0.0 {
            panic!("b/T must be non-negative. {} is invalid", b_by_tt);
        }
        let x = 0.5 * PI * s_by_tt;
        let y = 0.5 * PI * b_by_tt;
        let beta = f64::cos(x) * f64::sqrt(f64::powi(f64::tanh(y), 2) + f64::powi(f64::tan(x), 2));
        let beta_dash = f64::sqrt(1.0 - f64::powi(beta, 2));
        let kk_beta = elliptic_f(COMPLETE, beta).unwrap();
        let kk_beta_dash = elliptic_f(COMPLETE, beta_dash).unwrap();
        0.5 * kk_beta_dash / kk_beta
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod test {
    use super::ConfinedFlow;
    use plotpy::{linspace, Curve, Plot};
    use std::collections::HashMap;

    const SAVE_FIGURE: bool = true;

    #[test]
    fn test_confined_flow() {
        // let m = 0.655794203;
        // let kk_approx = |x: f64| { 0.512 / f64::sqrt(1.0 - 0.776 * x) + 0.258 / f64::sqrt(1.0 - 0.987 * x) + 0.8 / f64::sqrt(1.0 - 0.177 * x) };
        // let kk_beta = elliptic_f(PI / 2.0, m).unwrap();
        // let kk_beta_approx = kk_approx(m);
        // println!( "{} =? {} ({})", kk_beta, kk_beta_approx, f64::abs(kk_beta - kk_beta_approx));

        // NOTE: the reference data comes from a digitized plot probably from a scanned figure.
        // So, the results from the digitized plot may not be very accurate.

        // Reference data for b/T = 0.00 (x = s/T, y = Q/(kh))
        let ref_000_x = [
            0.01881, 0.02436, 0.03462, 0.04488, 0.06077, 0.07386, 0.09541, 0.11977, 0.14133, 0.17135, 0.20702, 0.24363,
            0.28026, 0.32441, 0.36576, 0.41839, 0.48418, 0.56501, 0.64961, 0.72294, 0.78404, 0.82164, 0.85266, 0.87616,
            0.90436, 0.92785, 0.94852, 0.96542, 0.97951, 0.98513, 0.99546,
        ];
        let ref_000_y = [
            1.4925, 1.42308, 1.34897, 1.27861, 1.20075, 1.13508, 1.06848, 0.9878, 0.92683, 0.8546, 0.79362, 0.73264,
            0.68292, 0.6379, 0.60413, 0.56004, 0.51407, 0.45685, 0.39962, 0.36022, 0.32364, 0.30206, 0.28518, 0.26923,
            0.24953, 0.22701, 0.20356, 0.18292, 0.16135, 0.1454, 0.12476,
        ];

        // Reference data for b/T = 0.25 (x = s/T, y = Q/(kh))
        let ref_025_x = [
            0.00188, 0.03574, 0.07148, 0.10251, 0.12883, 0.15798, 0.18429, 0.21061, 0.23129, 0.25385, 0.27452, 0.29332,
            0.31494, 0.33844, 0.36288, 0.38356, 0.40518, 0.42868, 0.45312, 0.47568, 0.49824, 0.51986, 0.54336, 0.56498,
            0.59036, 0.6195, 0.64394, 0.66744, 0.69095, 0.71351, 0.73513, 0.76145, 0.78401, 0.80751, 0.83102, 0.85452,
            0.87896, 0.90152, 0.92032, 0.93911, 0.95414, 0.96917, 0.98231, 0.99169,
        ];
        let ref_025_y = [
            0.75234, 0.74484, 0.73546, 0.72233, 0.70919, 0.69512, 0.67636, 0.66135, 0.6454, 0.62852, 0.61444, 0.5985,
            0.58536, 0.56942, 0.55535, 0.54315, 0.52908, 0.51313, 0.49906, 0.48499, 0.47185, 0.45685, 0.44277, 0.43058,
            0.41369, 0.39681, 0.3818, 0.36866, 0.35647, 0.3424, 0.32927, 0.31613, 0.30206, 0.28893, 0.2758, 0.26079,
            0.24765, 0.23546, 0.22138, 0.20263, 0.18855, 0.17073, 0.14728, 0.12195,
        ];

        // Reference data for b/T = 0.50 (x = s/T, y = Q/(kh))
        let ref_050_x = [
            0.00068, 0.03642, 0.07216, 0.10508, 0.14082, 0.17092, 0.20383, 0.23863, 0.27813, 0.31198, 0.34395, 0.37875,
            0.41072, 0.44081, 0.47654, 0.51039, 0.54612, 0.57903, 0.61382, 0.64484, 0.67869, 0.7116, 0.74544, 0.77553,
            0.80467, 0.82536, 0.84698, 0.86765, 0.89492, 0.91559, 0.93345, 0.95036, 0.96351, 0.98135,
        ];
        let ref_050_y = [
            0.5394, 0.53564, 0.53096, 0.52533, 0.51876, 0.51313, 0.50563, 0.49625, 0.48499, 0.47467, 0.46247, 0.45403,
            0.44184, 0.43152, 0.41932, 0.40525, 0.39118, 0.37711, 0.36304, 0.34896, 0.3349, 0.31895, 0.303, 0.2908,
            0.27485, 0.2636, 0.25047, 0.23452, 0.22045, 0.20263, 0.18762, 0.17261, 0.15572, 0.12664,
        ];

        // Reference data for b/T = 0.75 (x = s/T, y = Q/(kh))
        let ref_075_x = [
            0.00053, 0.05603, 0.10777, 0.16421, 0.23286, 0.30529, 0.38711, 0.46047, 0.53476, 0.60529, 0.66454, 0.73224,
            0.79524, 0.84507, 0.87703, 0.91181, 0.93719, 0.95127, 0.97006,
        ];
        let ref_075_y = [
            0.42495, 0.42307, 0.41838, 0.41369, 0.39962, 0.38837, 0.37148, 0.35553, 0.33583, 0.31707, 0.30112, 0.28049,
            0.25609, 0.23077, 0.21013, 0.18855, 0.16885, 0.14728, 0.12476,
        ];

        // Reference data for b/T = 1.00 (x = s/T, y = Q/(kh))
        let ref_100_x = [
            0.00045, 0.05783, 0.10486, 0.166, 0.22432, 0.28358, 0.35694, 0.42372, 0.49145, 0.55634, 0.62312, 0.68801,
            0.74256, 0.79239, 0.83846, 0.87701, 0.91178, 0.94186, 0.96159, 0.97378,
        ];
        let ref_100_y = [
            0.35835, 0.35366, 0.3499, 0.34521, 0.33771, 0.33396, 0.32364, 0.31238, 0.30488, 0.29268, 0.28049, 0.26735,
            0.2514, 0.22983, 0.21013, 0.19043, 0.16791, 0.1454, 0.12195, 0.09568,
        ];

        // Reference data for b/T = 1.25 (x = s/T, y = Q/(kh))
        let ref_125_x = [
            0.00226, 0.06153, 0.11044, 0.17535, 0.24872, 0.30421, 0.37664, 0.46787, 0.54406, 0.61084, 0.6701, 0.73594,
            0.78484, 0.84408, 0.89391, 0.92869, 0.95123, 0.97377,
        ];
        let ref_125_y = [
            0.30488, 0.30112, 0.29831, 0.29549, 0.28893, 0.28424, 0.27298, 0.26079, 0.25422, 0.2439, 0.23639, 0.22608,
            0.21294, 0.18761, 0.16698, 0.14634, 0.11913, 0.08724,
        ];

        // Reference data for b/T = 1.50 (x = s/T, y = Q/(kh))
        let ref_150_x = [
            0.00221, 0.05489, 0.10851, 0.17248, 0.22327, 0.28724, 0.36532, 0.44433, 0.5224, 0.59859, 0.66255, 0.73873,
            0.80174, 0.85816, 0.90611, 0.93054, 0.95779, 0.97564, 0.98502, 0.9944, 1.0,
        ];
        let ref_150_y = [
            0.25985, 0.25797, 0.25797, 0.25609, 0.25516, 0.25234, 0.24953, 0.2439, 0.23639, 0.22608, 0.2167, 0.20262,
            0.1848, 0.1651, 0.14071, 0.12289, 0.0985, 0.07692, 0.05534, 0.03283, 0.0075,
        ];

        let reference = HashMap::from([
            ("000", (ref_000_x.as_ref(), ref_000_y.as_ref())),
            ("025", (ref_025_x.as_ref(), ref_025_y.as_ref())),
            ("050", (ref_050_x.as_ref(), ref_050_y.as_ref())),
            ("075", (ref_075_x.as_ref(), ref_075_y.as_ref())),
            ("100", (ref_100_x.as_ref(), ref_100_y.as_ref())),
            ("125", (ref_125_x.as_ref(), ref_125_y.as_ref())),
            ("150", (ref_150_x.as_ref(), ref_150_y.as_ref())),
        ]);
        let x = linspace(0.0, 1.0, 201);
        let mut y = vec![0.0; x.len()];
        let mut plot = Plot::new();
        for (b_by_tt, tol) in [
            (0.0, 0.17),
            (0.25, 0.03),
            (0.5, 0.03),
            (0.75, 0.03),
            (1.0, 0.05),
            (1.25, 0.05),
            (1.5, 0.07),
        ] {
            let key = format!("{:03}", (b_by_tt * 100.0));
            let (ref_x, ref_y) = reference.get(key.as_str()).unwrap();
            for (i, s_by_tt) in ref_x.iter().enumerate() {
                let q_by_kh = ConfinedFlow::polubarinova_kochina_solution(*s_by_tt, b_by_tt);
                let q_ref = ref_y[i];
                let diff = f64::abs(q_by_kh - q_ref);
                assert!(
                    diff < tol,
                    "b/T={} s/T={} computed Q/(kH)={} reference Q/(kH)={} diff={}",
                    b_by_tt,
                    s_by_tt,
                    q_by_kh,
                    q_ref,
                    diff
                );
            }
            if SAVE_FIGURE {
                let mut curve = Curve::new();
                let mut curve_ref = Curve::new();
                curve_ref.set_marker_style("+").set_line_style("None");
                for i in 0..x.len() {
                    println!(
                        "b/T={}, s/T={} computed Q/(kH)={}",
                        b_by_tt,
                        x[i],
                        ConfinedFlow::polubarinova_kochina_solution(x[i], b_by_tt)
                    );
                    y[i] = ConfinedFlow::polubarinova_kochina_solution(x[i], b_by_tt);
                }
                curve.set_label(&format!("b/T={:.2}", b_by_tt));
                curve_ref.draw(ref_x, ref_y);
                if b_by_tt == 0.0 {
                    curve.draw(&x[1..].to_vec(), &y[1..].to_vec());
                } else {
                    curve.draw(&x, &y);
                }
                plot.add(&curve).add(&curve_ref);
            }
        }
        if SAVE_FIGURE {
            plot.set_title("Confined Seepage Flow - Polubarinova-Kochina Solution")
                .grid_labels_legend("s/T", "q/(kH)")
                .set_ymax(1.5)
                .set_figure_size_points(600.0, 800.0)
                .save("/tmp/pmsim/confined_seepage_flow_polubarinova_kochina_solution.svg")
                .unwrap();
        }
    }
}
