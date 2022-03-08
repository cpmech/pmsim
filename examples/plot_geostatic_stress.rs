use gemlab::mesh::Mesh;
use plotpy::{Curve, Plot, Shapes};
use pmsim::*;
use russell_lab::Vector;
use std::path::Path;

const OUT_DIR: &str = "/tmp/pmsim/examples";

fn main() -> Result<(), StrError> {
    let k_iso = 1e-2; // m/s
    let kk0 = 0.5;
    let nu = kk0 / (kk0 + 1.0);
    println!("nu = {}", nu);

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
        porosity_initial: 0.4,
        density_solid: 2.7, // Mg/m³
        stress_strain,
        retention_liquid,
        conductivity_liquid,
        conductivity_gas: None,
    };

    let lower = ParamPorous {
        earth_pres_coef_ini: kk0,
        porosity_initial: 0.1,
        density_solid: 2.7, // Mg/m³
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

    let mesh = Mesh::from_text_file("./data/meshes/column_two_layers_quads.msh")?;

    let mut config = SimConfig::new(&mesh);
    config
        .elements(1, ElementConfig::Porous(lower, None))?
        .elements(2, ElementConfig::Porous(upper, None))?
        .set_param_fluids(fluids)?
        .set_gravity(10.0)?; // m/s²

    let geo = Geostatics::new(&config)?;

    let count = 13;
    let zz = Vector::linspace(0.0, 3.0, count)?;
    let mut pl = Vector::new(count);
    let mut sig_h_eff = Vector::new(count);
    let mut sig_v_eff = Vector::new(count);
    let mut sig_h_tot = Vector::new(count);
    let mut sig_v_tot = Vector::new(count);

    for i in 0..count {
        let sig_eff = geo.calc_stress(zz[i], false)?;
        let sig_tot = geo.calc_stress(zz[i], true)?;
        pl[i] = geo.calc_pl(zz[i])?;
        sig_h_eff[i] = -sig_eff.get(0, 0); // negative to convert to soil mechanics' convention
        sig_v_eff[i] = -sig_eff.get(1, 1);
        sig_h_tot[i] = -sig_tot.get(0, 0);
        sig_v_tot[i] = -sig_tot.get(1, 1);
    }

    let mut curve_pl = Curve::new();
    let mut curve_sig_h_eff = Curve::new();
    let mut curve_sig_v_eff = Curve::new();
    let mut curve_sig_h_tot = Curve::new();
    let mut curve_sig_v_tot = Curve::new();
    curve_pl.set_label("pl").draw(&pl, &zz);
    curve_sig_h_eff
        .set_line_style("--")
        .set_label("sig_h_eff")
        .draw(&sig_h_eff, &zz);
    curve_sig_v_eff
        .set_line_style("--")
        .set_label("sig_v_eff")
        .draw(&sig_v_eff, &zz);
    curve_sig_h_tot
        .set_label("sig_h_tot")
        .set_line_width(2.0)
        .draw(&sig_h_tot, &zz);
    curve_sig_v_tot
        .set_label("sig_v_tot")
        .set_line_width(2.0)
        .draw(&sig_v_tot, &zz);

    let mut hl = Shapes::new();
    hl.set_edge_color("black")
        .draw_polyline(&vec![vec![0.0, 1.0], vec![70.0, 1.0]], false);

    let mut plot = Plot::new();
    plot.add(&hl)
        .add(&curve_pl)
        .add(&curve_sig_h_eff)
        .add(&curve_sig_v_eff)
        .add(&curve_sig_h_tot)
        .add(&curve_sig_v_tot)
        .grid_and_labels("pressure or stress component", "elevation")
        .set_figure_size_points(600.0, 600.0)
        .legend();

    let path = Path::new(OUT_DIR).join("plot_geostatic_stress.svg");
    plot.save(&path)?;
    Ok(())
}
