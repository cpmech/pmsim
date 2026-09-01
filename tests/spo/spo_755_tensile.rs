use gemlab::prelude::*;
use plotpy::{Curve, DarkMode, Legend, Plot, Text};
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;
use russell_lab::math::{PI, SQRT_3};
use russell_lab::{approx_eq, read_data, Norm, Vector};

const DIR: &str = "/tmp/pmsim/spo/spo_755";
const MESH_NAME: &str = "spo_755_tensile";
const NAME: &str = "spo_755_tensile";
const DRAW_MESH_AND_EXIT: bool = false;
const VERBOSE_LEVEL: usize = 0;
const SAVE_FIGURE: bool = false;

const YOUNG: f64 = 206.9; // Young's modulus
const POISSON: f64 = 0.29; // Poisson's coefficient
const KAPPA_INI: f64 = 0.45; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points
const WIDTH: f64 = 10.0; // width of the specimen
const B: f64 = 1.0; // width of the ligament

// loading factors
const LAMBDAS: [f64; 17] = [
    0.0,   //  0
    0.005, //  1
    0.01,  //  2
    0.015, //  3
    0.02,  //  4
    0.03,  //  5
    0.04,  //  6
    0.05,  //  7
    0.07,  //  8
    0.09,  //  9
    0.10,  //  9b (need this extra step compared to SPO's code)
    0.11,  // 10
    0.115, // 11
    0.12,  // 12
    0.13,  // 13
    0.15,  // 14
    0.17,  // 15
];

#[test]
fn spo_755_tensile() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", MESH_NAME))?;

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let top = features.search_edges(At::Y(15.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), |x| x[0] < 0.50001)?;

    // draw mesh
    if DRAW_MESH_AND_EXIT {
        let ids_left = features.get_points_via_2d_edges(&left);
        let ids_top = features.get_points_via_2d_edges(&top);
        let ids_bottom = features.get_points_via_2d_edges(&bottom);
        println!("ids left = {:?}", ids_left);
        println!("ids top = {:?}", ids_top);
        println!("ids bottom = {:?}", ids_bottom);
        return draw_mesh(&mesh, &ids_left, &ids_bottom, &ids_top);
    }

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            kappa_ini: KAPPA_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };

    // schema
    let mut schema = Schema::new();
    schema.add_solid(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges(&top, Dof::Uy, 1.0);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // run: natural + sps
    let options = Options {
        arclength: false,
        lmm: false,
    };
    run(options, &mesh, &features, &bottom, &schema, &mut ebc, &nbc)?;
    Ok(())
}

// simulation ----------------------------------------------------------------
fn run(
    options: Options,
    mesh: &Mesh,
    features: &Features,
    bottom: &Edges,
    schema: &Schema,
    ebc: &mut BcEssential,
    nbc: &BcNatural,
) -> Result<(), StrError> {
    // define filename stem
    let mut name = NAME.to_string() + "_";
    name += &options.key();

    // find corner node
    let (min, max) = mesh.get_limits();
    let corner_id = features.search_point_ids(At::XY(min[0], max[1]), any_x)?[0];

    // configuration
    let mut config = Config::<2>::new(&mesh);
    config
        .out_history_uu_comp(corner_id, Dof::Uy)
        .lagrange_mult_method(options.lmm);

    // output files if natural parameter continuation (for verification)
    if !options.arclength {
        config.out_files(DIR, &name);
    }

    // output the vertical component of Y at bottom edge points
    let ids_bottom = features.get_points_via_2d_edges(&bottom);
    for point_id in &ids_bottom {
        config.out_history_yy_comp(*point_id, Dof::Uy);
    }

    // nonlinear solver configuration
    let mut nl_config = NlConfig::new();
    nl_config
        .set_verbose(false, true, true)
        .set_record_iterations_residuals(false)
        .set_tg_control_tol(0.5);
    if options.arclength {
        nl_config.set_method(NlMethod::Arclength).set_bordering(true);
    }

    // simulator and data
    let (mut sim, mut data) = Simulator::new(&mesh, &schema, &config, &ebc, &nbc, &mut nl_config)?;

    // stopping criteria
    let stop = if options.arclength {
        // note: cannot use the corner point on the SPS because only unknown equations are available
        let end = usize::min(data.nsys(), data.ndof()); // skip Lagrange multipliers, if any
        Stop::MaxNormU(0.25, Norm::Max, 0, end)
    } else {
        Stop::Steps(LAMBDAS.len() - 1)
    };

    // delta lambda
    let dll = if options.arclength {
        DeltaLambda::auto(0.005)
    } else {
        let list = Vector::from(&LAMBDAS).get_differences();
        DeltaLambda::list(list.as_data())
    };

    // run simulation
    sim.steady(&mut data, IniDir::Pos, stop, dll)?;

    //
    // data analysis -------------------------------------------------------------
    //

    // calculate the reaction using the internal forces
    let lambdas = data.stations();
    if !options.arclength {
        assert_eq!(lambdas, &LAMBDAS);
    }
    let nstation = lambdas.len();
    let mut sum_yy = vec![0.0; nstation];
    for point_id in &ids_bottom {
        let yy_over_time = data.history_yy_comp(*point_id, Dof::Uy).unwrap();
        for i in 0..nstation {
            sum_yy[i] += yy_over_time[i];
        }
    }

    // get the history of vertical displacement at the corner point
    let history_uy = data.history_uu_comp(corner_id, Dof::Uy).unwrap();

    // loop over stations (lambdas)
    let analytical_limit = (2.0 + PI) / SQRT_3;
    let mut normalized_deflection = Vec::with_capacity(nstation);
    let mut normalized_stress = Vec::with_capacity(nstation);
    for index in 0..nstation {
        let uy = history_uy[index];
        normalized_deflection.push(2.0 * uy * YOUNG / (KAPPA_INI * WIDTH));
        let reaction = 2.0 * sum_yy[index]; // multiply by 2 because only half specimen is modeled
        let tensile_stress = -reaction / B; // negative because the reaction points downwards
        let norm_net_stress = tensile_stress / KAPPA_INI;
        normalized_stress.push(norm_net_stress);
        if index == nstation - 1 {
            let diff = f64::abs(norm_net_stress - analytical_limit);
            println!(
                "final normalized net stress = {}, diff = {} ({:.2}%)",
                norm_net_stress,
                diff,
                100.0 * diff / analytical_limit
            );
            let tol = if options.arclength {
                if options.lmm {
                    0.0255
                } else {
                    0.027
                }
            } else {
                0.02535
            };
            approx_eq(norm_net_stress, analytical_limit, tol);
        }
    }

    // plot the results
    if SAVE_FIGURE {
        let title = options.title();
        let ref_xy = read_data("data/spo/spo_755_fig730a.tsv", &["x", "y"])?;
        let mut plot = Plot::new();
        let mut curve_num = Curve::new();
        let mut curve_ref = Curve::new();
        curve_ref
            .set_label("de Souza Neto et al. (scanned)")
            .set_line_style(":")
            .set_line_color("#1ea56a")
            .set_marker_style("+")
            .set_marker_size(12.0)
            .draw(&ref_xy["x"], &ref_xy["y"]);
        curve_num
            .set_label("pmsim")
            .set_line_color("#8e0220")
            .set_marker_style("o")
            .draw(&normalized_deflection, &normalized_stress);
        let mut txt = Text::new();
        txt.set_align_horizontal("left").set_align_vertical("bottom").draw(
            0.0,
            analytical_limit,
            "$(2 + \\pi)/\\sqrt{3}$",
        );
        let mut dm = DarkMode::new();
        let mut leg = Legend::new();
        dm.set_mocha();
        leg.set_location("lower right").draw();
        plot.add(&dm)
            .set_horiz_line(analytical_limit, "#51b4df", "--", 1.0)
            .add(&curve_ref)
            .add(&curve_num)
            .add(&leg)
            .add(&txt)
            .set_rotation_ticks_x(90.0)
            .grid_and_labels(
                "$2 u_y E / (\\sigma_y w)$ (normalized deflection)",
                "$\\bar{\\sigma}/\\sigma_y$ (normalized stress)",
            )
            .set_figure_size_points(600.0, 600.0)
            .set_title(&title)
            .save(&format!("{}/{}.svg", DIR, &name))
            .unwrap();
    }

    //
    // verification --------------------------------------------------------------
    //

    // verify the results
    if !options.arclength {
        let tol_displacement = 3.38e-9;
        let tol_stress = 1.13e-5;
        let all_good = compare_results(
            &mesh,
            &schema,
            &config,
            DIR,
            &name,
            ReferenceDataType::SPO,
            &format!("data/spo/{}_ref.json", NAME),
            tol_displacement,
            tol_stress,
            VERBOSE_LEVEL,
            Some((1, 1.0, 1e-6)),
        )?;
        assert!(all_good);
    }
    Ok(())
}

fn draw_mesh(mesh: &Mesh, left: &[PointId], bottom: &[PointId], top: &[PointId]) -> Result<(), StrError> {
    mesh.check_all()?;
    let spo_fixities = [
        (1, "10"),
        (2, "10"),
        (5, "01"),
        (6, "11"),
        (7, "10"),
        (9, "01"),
        (15, "01"),
        (19, "01"),
        (30, "11"),
        (33, "10"),
        (36, "01"),
        (37, "10"),
        (38, "10"),
        (53, "01"),
        (54, "01"),
        (63, "01"),
        (64, "01"),
        (81, "10"),
        (82, "10"),
        (92, "10"),
        (93, "10"),
        (181, "10"),
        (184, "10"),
        (222, "01"),
        (224, "10"),
        (225, "01"),
        (230, "10"),
        (235, "10"),
        (253, "10"),
        (260, "10"),
        (266, "10"),
        (269, "01"),
        (279, "01"),
        (286, "01"),
        (295, "01"),
        (302, "01"),
        (366, "10"),
        (373, "10"),
        (380, "10"),
        (449, "01"),
        (451, "01"),
        (31, "01"),
        (32, "01"),
        (83, "01"),
        (84, "01"),
        (85, "01"),
        (86, "01"),
        (379, "01"),
        (382, "01"),
        (384, "01"),
        (469, "01"),
        (471, "01"),
        (473, "01"),
    ];
    let mut text = Text::new();
    text.set_extra("clip_on=True")
        .set_fontsize(6.0)
        .set_align_horizontal("center")
        .set_align_vertical("bottom")
        .set_bbox(true)
        .set_bbox_facecolor("yellow")
        .set_bbox_edgecolor("None")
        .set_bbox_style("round,pad=0.1,rounding_size=0.15");
    let mut spo_left = Vec::new();
    let mut spo_top = Vec::new();
    let mut spo_bottom = Vec::new();
    for (spo_id, fix) in spo_fixities.iter() {
        let id = spo_id - 1;
        let p = &mesh.points[id];
        text.draw(p.coords[0], p.coords[1], &format!("{}", fix));
        if f64::abs(p.coords[0]) < 0.0001 {
            spo_left.push(id);
        }
        if f64::abs(p.coords[1]) > 14.9999 {
            spo_top.push(id);
        }
        if f64::abs(p.coords[1]) < 0.0001 {
            spo_bottom.push(id);
        }
    }
    spo_left.sort();
    spo_top.sort();
    spo_bottom.sort();
    assert_eq!(&left, &spo_left);
    assert_eq!(&top, &spo_top);
    assert_eq!(&bottom, &spo_bottom);
    let mut draw = Draw::new();
    draw.set_range_2d(-0.5, 15.5, -0.5, 15.5)
        .set_size(1200.0, 1200.0)
        .zoom_2d(-0.02, 0.52, -0.02, 0.52, 0.38, 0.1, 0.6, 0.6)
        .zoom_extra(|inset| {
            inset.add(&text);
        })
        .extra(|plot, before| {
            if !before {
                plot.add(&text);
            }
        })
        .all(&mesh, &format!("{}/{}_mesh.svg", DIR, NAME))
}

struct Options {
    arclength: bool,
    lmm: bool,
}

impl Options {
    fn key(&self) -> String {
        let mut buf = if self.arclength {
            "arc".to_string()
        } else {
            "nat".to_string()
        };
        if self.lmm {
            buf += "_lmm";
        } else {
            buf += "_sps";
        }
        buf
    }

    fn title(&self) -> String {
        self.key().to_uppercase().replace("_", " | ")
    }
}
