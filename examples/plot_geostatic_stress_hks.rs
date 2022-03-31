use gemlab::mesh::Mesh;
use plotpy::{Curve, Plot, Shapes, Text};
use pmsim::{geostatics::*, simulation::*, StrError};
use russell_chk::assert_approx_eq;
use russell_lab::Vector;
use std::path::Path;

const OUT_DIR: &str = "/tmp/pmsim/examples";

fn print_values_and_check(geo: &Geostatics) -> Result<(), StrError> {
    let sig_t_a = geo.calc_sigma_z_total(6.0)?;
    let sig_t_b = geo.calc_sigma_z_total(4.0001)?;
    let sig_t_c = geo.calc_sigma_z_total(4.0)?;
    let sig_t_d = geo.calc_sigma_z_total(3.9999)?;
    let sig_t_e = geo.calc_sigma_z_total(0.0)?;
    let pl_a = geo.calc_pl(6.0)?;
    let pl_b = geo.calc_pl(4.0001)?;
    let pl_c = geo.calc_pl(4.0)?;
    let pl_d = geo.calc_pl(3.9999)?;
    let pl_e = geo.calc_pl(0.0)?;
    let sig_e_a = geo.calc_stress(6.0, false)?.vec[1];
    let sig_e_b = geo.calc_stress(4.0001, false)?.vec[1];
    let sig_e_c = geo.calc_stress(4.0, false)?.vec[1];
    let sig_e_d = geo.calc_stress(3.9999, false)?.vec[1];
    let sig_e_e = geo.calc_stress(0.0, false)?.vec[1];
    println!("            z        sig_v    pl  sig_v'");
    println!("surface:    6.0    {:>7.2} {:>5.2} {:>5.2}", sig_t_a, pl_a, sig_e_a);
    println!("transition: 4.0001 {:>7.2} {:>5.2} {:>5.2}", sig_t_b, pl_b, sig_e_b);
    println!("transition: 4.0    {:>7.2} {:>5.2} {:>5.2}", sig_t_c, pl_c, sig_e_c);
    println!("transition: 3.9999 {:>7.2} {:>5.2} {:>5.2}", sig_t_d, pl_d, sig_e_d);
    println!("bottom:     0.0    {:>7.2} {:>5.2} {:>5.2}", sig_t_e, pl_e, sig_e_e);
    assert_approx_eq!(sig_t_a, -26.49, 1e-15);
    assert_approx_eq!(sig_t_b, -62.79, 1e-2);
    assert_approx_eq!(sig_t_c, -62.79, 1e-2);
    assert_approx_eq!(sig_t_d, -62.79, 1e-2);
    assert_approx_eq!(sig_t_e, -141.27, 1e-2);
    assert_approx_eq!(pl_a, 0.0, 1e-15);
    assert_approx_eq!(pl_b, 19.62, 1e-3);
    assert_approx_eq!(pl_c, 19.62, 1e-3);
    assert_approx_eq!(pl_d, 19.62, 1e-3);
    assert_approx_eq!(pl_e, 58.86, 1e-3);
    assert_approx_eq!(sig_e_a, -26.49, 1e-15);
    assert_approx_eq!(sig_e_b, -43.17, 1e-2);
    assert_approx_eq!(sig_e_c, -43.17, 1e-2);
    assert_approx_eq!(sig_e_d, -43.17, 1e-2);
    assert_approx_eq!(sig_e_e, -82.41, 1e-2);
    Ok(())
}

fn main() -> Result<(), StrError> {
    // Examples 6.8 and 6.9 (pages 271-273) of the following book:
    //
    // * Robert D. Holtz, William D. Kovacs, and Thomas C. Sheehan,
    //   An Introduction to Geotechnical Engineering, Pearson
    //
    // The top-most layer (unsaturated sand) in the book is modelled by an
    // equivalent overburden stress, thus, here we only consider two layers.
    //
    // The overburden given by the unsaturated sand layer is equal to -26.49 kN/m²
    // (negative means compression)
    //
    // Assuming nf=0.4 for the clay layer (lower), thus:
    //
    //      ρSat - ρW・nf   2 - 1・0.4
    // ρS = ————————————— = ————————— = 2.666666666666667
    //         1 - nf        1 - 0.4

    let k_iso = 1e-2; // m/s
    let kk0 = 0.25; // assumed
    let nu = kk0 / (kk0 + 1.0);

    let stress_strain = ParamStressStrain::LinearElastic {
        young: 10_000.0, // kPa
        poisson: nu,     // [-]
    };

    let retention_liquid = ParamLiquidRetention::PedrosoWilliams {
        with_hysteresis: true,
        lambda_d: 3.0,
        lambda_w: 3.0,
        beta_d: 6.0,
        beta_w: 6.0,
        beta_1: 6.0,
        beta_2: 6.0,
        x_rd: 2.0,
        x_rw: 2.0,
        y_0: 1.0,
        y_r: 0.005,
    };

    let conductivity_liquid = ParamConductivity::PedrosoZhangEhlers {
        kx: k_iso, // m/s
        ky: k_iso, // m/s
        kz: k_iso, // m/s
        lambda_0: 0.001,
        lambda_1: 1.2,
        alpha: 0.01,
        beta: 10.0,
    };

    let upper = ParamPorous {
        earth_pres_coef_ini: kk0,
        porosity_initial: 0.5,
        density_solid: 2.7, // Mg/m³
        stress_strain,
        retention_liquid,
        conductivity_liquid,
        conductivity_gas: None,
    };

    let lower = ParamPorous {
        earth_pres_coef_ini: kk0,
        porosity_initial: 0.4,
        density_solid: 2.666666666666667, // Mg/m³
        stress_strain,
        retention_liquid,
        conductivity_liquid,
        conductivity_gas: None,
    };

    let fluids = ParamFluids {
        density_liquid: ParamRealDensity {
            cc: 1e-12,    // Mg/(m³ kPa)
            p_ref: 0.0,   // kPa
            rho_ref: 1.0, // Mg/m³
            tt_ref: 25.0, // ℃
        },
        density_gas: None,
    };

    let mesh = Mesh::from_text_file("./data/meshes/column_hks_example.msh")?;

    let mut config = Configuration::new(&mesh);
    config
        .elements(1, ElementConfig::Porous(lower, None))?
        .elements(2, ElementConfig::Porous(upper, None))?
        .fluids(fluids)?
        .gravity(9.81)? // m/s²
        .ini_option(IniOption::Geostatic(-26.49))?; // kN/m²

    let geo = Geostatics::new(&config)?;

    print_values_and_check(&geo)?;

    let count = 41;
    let zz = Vector::linspace(0.0, 6.0, count)?;
    let mut pl = Vector::new(count);
    let mut sig_v_eff = Vector::new(count);
    let mut sig_v_tot = Vector::new(count);

    for i in 0..count {
        let sig_eff = geo.calc_stress(zz[i], false)?;
        let sig_tot = geo.calc_stress(zz[i], true)?;
        pl[i] = geo.calc_pl(zz[i])?;
        sig_v_eff[i] = -sig_eff.get(1, 1); // negative to convert to soil mechanics' convention
        sig_v_tot[i] = -sig_tot.get(1, 1);
    }

    let mut curve_pl = Curve::new();
    let mut curve_sig_v_eff = Curve::new();
    let mut curve_sig_v_tot = Curve::new();
    curve_pl.set_label("pl").draw(&pl, &zz);
    curve_sig_v_eff
        .set_line_style("--")
        .set_label("sig_v_eff")
        .draw(&sig_v_eff, &zz);
    curve_sig_v_tot
        .set_label("sig_v_tot")
        .set_line_width(2.0)
        .draw(&sig_v_tot, &zz);

    let mut hl = Shapes::new();
    hl.set_edge_color("black")
        .draw_polyline(&vec![vec![0.0, 4.0], vec![150.0, 4.0]], false);

    let mut t1 = Text::new();
    let mut t2 = Text::new();
    let mut t3 = Text::new();
    let mut t4 = Text::new();
    let mut t5 = Text::new();
    let mut t6 = Text::new();
    let mut t7 = Text::new();
    t1.draw(26.49, 6.0, "26.49");
    t2.draw(62.79, 4.0, "62.79");
    t3.draw(141.27, 0.0, "141.27");
    t4.draw(43.17, 4.0, "43.17");
    t5.draw(82.41, 0.0, "82.41");
    t6.draw(19.62, 4.0, "19.62");
    t7.draw(58.86, 0.0, "58.86");

    let mut plot = Plot::new();
    plot.add(&hl)
        .add(&curve_pl)
        .add(&curve_sig_v_eff)
        .add(&curve_sig_v_tot)
        .add(&t1)
        .add(&t2)
        .add(&t3)
        .add(&t4)
        .add(&t5)
        .add(&t6)
        .add(&t7)
        .grid_and_labels("pressure or stress component", "elevation")
        .set_figure_size_points(600.0, 600.0)
        .legend();

    let path = Path::new(OUT_DIR).join("plot_geostatic_stress_hks.svg");
    plot.save(&path)?;
    Ok(())
}
