use gemlab::prelude::*;
use pmsim::fem::PostProcMemo;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{vec_inner, Vector};

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_cofferdam";
const SAVE_FIGURE: bool = true;

#[test]
fn test_seep_cofferdam() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read("data/meshes/cofferdam.msh")?;
    if SAVE_FIGURE {
        let x_ticks = [-13.6, -4.6, 0.0, 4.4, 13.6];
        let x_labels = x_ticks.iter().map(|x| format!("{}", x)).collect::<Vec<_>>();
        let y_ticks = [0.0, 5.0, 7.5, 10.0];
        let y_labels = y_ticks.iter().map(|y| format!("{}", y)).collect::<Vec<_>>();
        let mut draw = Draw::new();
        draw.show_edge_markers(true)
            .show_point_ids(true)
            .show_cell_ids(true)
            .set_m_normal_vector_marker(0.025)
            .set_size(1200.0, 600.0)
            .set_view_flag(false)
            .extra(|plot, before| {
                if !before {
                    plot.set_ticks_x_labels(&x_ticks, &x_labels)
                        .set_ticks_y_labels(&y_ticks, &y_labels)
                        .set_range(-15.0, 15.0, -1.0, 11.0);
                }
            });
        draw.all(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))?;
    }

    // features
    let features = Features::new(&mesh, false);
    let inlet = features.search_marked_edges(20);
    let outlet = features.search_marked_edges(10);

    // parameters
    let (kx, ky) = (4e-7, 4e-7);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))])?;

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.edges(&inlet, Dof::Phi, 13.0);
    essential.edges(&outlet, Dof::Phi, 7.5);

    // natural boundary conditions
    let natural = Natural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_out_files(OUT_DIR, NAME, 1.0)
        .update_model_settings(1)
        .set_save_flux(true);

    // FEM state
    let mut state = FemState::new(&mesh, &base, &essential, &config)?;

    // FEM results
    let mut results = FemResults::new(&mesh, &base, &config)?;

    // solution
    let mut solver = SolverImplicit::new(&mesh, &base, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut results)?;

    // post-processing
    post_processing()
}

fn post_processing() -> Result<(), StrError> {
    // post-processing tool
    let (post, mut memo) = PostProc::new(OUT_DIR, NAME)?;

    // read last state
    let last = post.n_state() - 1;
    let state = post.read_state(last)?;

    // load mesh and find vertical section along the gap underneath the wall
    let mesh = post.mesh();
    let with_internal_edges = true; // need internal edges
    let cells_by_points = true; // use all cells surrounding a point for better extrapolation
    let features = Features::new(&mesh, with_internal_edges);
    let (_, top_cells, top_section) = features.search_features_2d(cells_by_points, At::Y(10.0), |x| x[0] <= -4.6)?;
    let (_, mid_cells, mid_section) = features.search_features_2d(cells_by_points, At::X(-4.6), |x| x[1] <= 5.0)?;
    let (_, rig_cells, rig_section) = features.search_features_2d(cells_by_points, At::X(4.6), |x| x[1] <= 5.0)?;
    if mesh.cells.len() < 30 {
        println!("top_cells = {:?}", top_cells);
        println!("mid_cells = {:?}", mid_cells);
        println!("rig_cells = {:?}", rig_cells);
        println!(
            "top_section = {:?}",
            top_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
        println!(
            "mid_section = {:?}",
            mid_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
        println!(
            "rig_section = {:?}",
            rig_section.all.iter().map(|e| e.key()).collect::<Vec<_>>()
        );
    }

    // flux through top section
    let un0 = Vector::from(&[0.0, -1.0]);
    let q0 = calc_flux_through_flat_section(&post, &mut memo, &state, &top_cells, &top_section, &un0)?;
    println!("Flux through top section = {:.3e}", q0);

    // flux through mid section
    let un1 = Vector::from(&[1.0, 0.0]);
    let q1 = calc_flux_through_flat_section(&post, &mut memo, &state, &mid_cells, &mid_section, &un1)?;
    println!("Flux through mid section = {:.3e}", q1);

    // flux through right section
    let un2 = Vector::from(&[-1.0, 0.0]); // to the left now
    let q2 = calc_flux_through_flat_section(&post, &mut memo, &state, &rig_cells, &rig_section, &un2)?;
    println!("Flux through rig section = {:.3e}", q2);

    println!("q1 + q2 = {:.3e}", q1 + q2);

    // write Paraview files
    let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    println!("Paraview File: {}", path_pvd);

    // done
    Ok(())
}

fn calc_flux_through_flat_section(
    post: &PostProc,
    memo: &mut PostProcMemo,
    state: &FemState,
    cells: &Vec<CellId>,
    section: &Edges,
    unit_normal: &Vector,
) -> Result<f64, StrError> {
    let patch = post.nodal_fluxes_patch(memo, &state, &cells, Dof::Phi, |_, _, _| true)?;
    let (points, coords, ww) = post.values_along_edges_vec(&patch, &section)?;
    let mut area = 0.0;
    let vertical_section = if unit_normal[0] == 0.0 {
        assert!(f64::abs(unit_normal[1]) == 1.0, "Unit normal must be unitary");
        false
    } else if unit_normal[1] == 0.0 {
        assert!(f64::abs(unit_normal[0]) == 1.0, "Unit normal must be unitary");
        true
    } else {
        panic!("Unit normal must be aligned with either x or y axis");
    };
    for i in 1..points.len() {
        let dx = coords[i][0] - coords[i - 1][0];
        let dy = coords[i][1] - coords[i - 1][1];
        if vertical_section {
            assert!(f64::abs(dx) < 1e-14, "Section must be vertical");
        } else {
            assert!(f64::abs(dy) < 1e-14, "Section must be horizontal");
        }
        let delta = if vertical_section { dy } else { dx };
        let w_dot_un_prev = vec_inner(&ww[i - 1], &unit_normal);
        let w_dot_un = vec_inner(&ww[i], &unit_normal);
        area += delta * (w_dot_un_prev + w_dot_un) / 2.0;
    }
    Ok(area)
}
