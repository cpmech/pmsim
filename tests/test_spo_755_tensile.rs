use gemlab::prelude::*;
use plotpy::Text;
use pmsim::prelude::*;
use pmsim::util::{compare_results, ReferenceDataType};
use pmsim::StrError;

const MESH_NAME: &str = "spo_755_tensile";
const NAME: &str = "spo_755_tensile_perf_plast";
const DRAW_MESH_AND_EXIT: bool = false;
const VERBOSE_LEVEL: usize = 0;

const YOUNG: f64 = 206.9; // Young's modulus
const POISSON: f64 = 0.29; // Poisson's coefficient
const Z_INI: f64 = 0.45; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_755_tensile() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/spo/{}.msh", MESH_NAME))?;

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let bottom = features.search_edges(At::Y(0.0), |x| x[0] < 0.50001)?;
    let top = features.search_edges(At::Y(15.0), any_x)?;

    // draw mesh
    if DRAW_MESH_AND_EXIT {
        println!("left = {:?}", features.get_points_via_2d_edges(&left));
        println!("bottom = {:?}", features.get_points_via_2d_edges(&bottom));
        println!("top = {:?}", features.get_points_via_2d_edges(&top));
        let ids_left = features.get_points_via_2d_edges(&left);
        let ids_bottom = features.get_points_via_2d_edges(&bottom);
        let ids_top = features.get_points_via_2d_edges(&top);
        return draw_mesh(&mesh, &ids_left, &ids_bottom, &ids_top);
    }

    // parameters
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: StressStrain::VonMises {
            young: YOUNG,
            poisson: POISSON,
            z_ini: Z_INI,
            hh: H,
        },
        ngauss: Some(NGAUSS),
    };
    let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))])?;

    // displacement control
    const UY_PERF_PLAST: [f64; 16] = [
        0.0, 0.005, 0.01, 0.015, 0.02, 0.03, 0.04, 0.05, 0.07, 0.09, 0.11, 0.115, 0.12, 0.13, 0.15, 0.17,
    ];

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&bottom, Dof::Uy, 0.0)
        .edges_fn(&top, Dof::Uy, 1.0, |t| UY_PERF_PLAST[t as usize]);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_incremental(UY_PERF_PLAST.len())
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, "/tmp/pmsim", NAME)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    // verify the results
    let tol_displacement = 1e-7;
    let tol_stress = 1.2e-5;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        ReferenceDataType::SPO,
        &format!("data/spo/{}_ref.json", NAME),
        tol_displacement,
        tol_stress,
        VERBOSE_LEVEL,
    )?;
    assert!(all_good);
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
    let mut spo_bottom = Vec::new();
    let mut spo_top = Vec::new();
    for (spo_id, fix) in spo_fixities.iter() {
        let id = spo_id - 1;
        let p = &mesh.points[id];
        text.draw(p.coords[0], p.coords[1], &format!("{}", fix));
        if f64::abs(p.coords[0]) < 0.0001 {
            spo_left.push(id);
        }
        if f64::abs(p.coords[1]) < 0.0001 {
            spo_bottom.push(id);
        }
        if f64::abs(p.coords[1]) > 14.9999 {
            spo_top.push(id);
        }
    }
    spo_left.sort();
    spo_bottom.sort();
    spo_top.sort();
    assert_eq!(&left, &spo_left);
    assert_eq!(&bottom, &spo_bottom);
    assert_eq!(&top, &spo_top);
    let mut fig = Figure::new();
    fig.range_2d(-0.5, 15.5, -0.5, 15.5)
        .size(1200.0, 1200.0)
        .zoom_2d(-0.02, 0.52, -0.02, 0.52, 0.38, 0.1, 0.6, 0.6)
        .zoom_extra(|inset| {
            inset.add(&text);
        })
        .extra(|plot, before| {
            if !before {
                plot.add(&text);
            }
        })
        .draw(&mesh, &format!("/tmp/pmsim/{}_mesh.svg", NAME))
}
