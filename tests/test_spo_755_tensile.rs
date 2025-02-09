use gemlab::prelude::*;
use plotpy::Text;
use pmsim::prelude::*;
use pmsim::util::compare_results;
use russell_lab::*;

const NAME: &str = "test_spo_755_tensile";
const DRAW_MESH: bool = true;

const YOUNG: f64 = 206.9; // Young's modulus
const POISSON: f64 = 0.29; // Poisson's coefficient
const Z_INI: f64 = 0.45; // Initial size of yield surface
const H: f64 = 0.0; // hardening coefficient
const NGAUSS: usize = 4; // number of gauss points

#[test]
fn test_spo_755_tensile() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read(&format!("data/meshes/{}.msh", NAME))?;
    if DRAW_MESH {
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
        ];
        let mut fig = Figure::new();
        fig.size(600.0, 1200.0)
            .extra(|plot, before| {
                if !before {
                    let mut text = Text::new();
                    text.set_bbox_style("round,pad=0.3,rounding_size=0.15");
                    for (spo_id, fix) in spo_fixities.iter() {
                        let id = spo_id - 1;
                        let p = &mesh.points[id];
                        text.draw(p.coords[0], p.coords[1], &format!("{}", fix));
                    }
                    plot.add(&text);
                }
            })
            .draw(&mesh, &format!("/tmp/pmsim/{}.svg", NAME))?;
        return Ok(());
    }

    // features
    let features = Features::new(&mesh, false);
    let left = features.search_edges(At::X(0.0), any_x)?;
    let notch = features.search_edges(At::Y(0.0), |x| x[0] >= 0.5)?;
    let top = features.search_edges(At::Y(15.0), any_x)?;
    println!("notch = {:?}", features.get_points_via_2d_edges(&notch));
    println!("top = {:?}", features.get_points_via_2d_edges(&top));

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
    const UY_PERF_PLAST: [f64; 5] = [
        0.0, 0.005, 0.01, 0.015, 0.02, // 0.03, 0.04, 0.05, 0.07, 0.09, 0.11, 0.115, 0.12, 0.13, 0.15, 0.17,
    ];

    // essential boundary conditions
    let mut essential = Essential::new();
    essential
        .edges(&left, Dof::Ux, 0.0)
        .edges(&notch, Dof::Uy, 0.0)
        .edges_fn(&top, Dof::Uy, |t| UY_PERF_PLAST[t as usize]);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_tol_rr(1e-7)
        .set_incremental(UY_PERF_PLAST.len())
        .set_symmetry_check_tolerance(Some(1e-5))
        .set_n_max_iterations(20);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // File IO
    let mut file_io = FileIo::new();
    file_io.activate(&mesh, &base, NAME, None)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut file_io)?;

    /*
    // verify the results
    let tol_displacement = 1e-15;
    let tol_stress = 1e-15;
    let all_good = compare_results(
        &mesh,
        &base,
        &file_io,
        "spo_755_tensile_perf_plast.json",
        tol_displacement,
        tol_stress,
        2,
    )?;
    assert!(all_good);
    */
    Ok(())
}
