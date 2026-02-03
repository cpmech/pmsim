use gemlab::prelude::*;
use pmsim::analytical::ConfinedFlow;
use pmsim::prelude::*;
use pmsim::StrError;
use russell_lab::{approx_eq, vec_inner, Vector};

const OUT_DIR: &str = "/tmp/pmsim";
const NAME: &str = "test_seep_simple_confined_flow_2d";
const SAVE_FIGURE: bool = false;

#[test]
fn test_seep_simple_confined_flow() -> Result<(), StrError> {
    // mesh
    let mesh = Mesh::read("data/meshes/simple-confined-flow-2d.msh")?;
    // let mesh = Mesh::read("data/meshes/simple-confined-flow-2d_fine.msh")?;
    if SAVE_FIGURE {
        let mut draw = Draw::new();
        draw.show_edge_markers(true)
            .show_point_ids(true)
            .set_m_normal_vector_marker(0.025)
            .set_size(1200.0, 600.0)
            .set_view_flag(false)
            .extra(|plot, before| {
                if !before {
                    plot.set_num_ticks_x(38)
                        .set_num_ticks_y(10)
                        .set_range(-1.0, 37.0, -1.0, 9.0);
                }
            });
        draw.all(&mesh, &format!("{}/mesh_{}.svg", OUT_DIR, NAME))?;
    }

    // features
    let with_internal_edges = true; // for data analysis
    let features = Features::new(&mesh, with_internal_edges);
    let inlet = features.search_marked_edges(12);
    let outlet = features.search_marked_edges(8);

    // parameters
    let (kx, ky) = (1.0, 1.0);
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: Conductivity::Constant { kx, ky, kz: 0.0 },
        source: None,
        ngauss: None,
    };
    let mut schema = Schema::new();
    schema.add_diffusion(1, p1).build(&mesh)?;

    // essential boundary conditions
    let mut ebc = BcEssential::new();
    ebc.edges(&inlet, Dof::Phi, 12.0);
    ebc.edges(&outlet, Dof::Phi, 8.0);

    // natural boundary conditions
    let nbc = BcNatural::new();

    // configuration
    let mut config = Config::new(&mesh);
    config
        .set_lagrange_mult_method(true)
        .set_out_files(OUT_DIR, NAME, 1.0)
        .update_model_settings(1)
        .set_save_flux(true);

    // solution
    let (mut sim, mut data) = SimulatorLin::new(&mesh, &schema, &config, &ebc, &nbc)?;
    sim.steady(&mut data, true)?;

    //
    // data analysis -------------------------------------------------------------
    //

    // post-processing tool
    let (post, mut memo) = PostProc::new(OUT_DIR, NAME)?;

    // load mesh and find vertical section along the gap underneath the wall
    let mid_section = features.search_edges(At::X(18.0), |_| true)?;

    // read last state
    let last = post.nstate() - 1;
    let state = post.read_state(last)?;

    // extract fluxes along the mid section
    let mid_cell_ids = features.get_cells_via_2d_edges(&mid_section);
    let mid_patch = post.nodal_fluxes_patch(&mut memo, &state, &mid_cell_ids, Dof::Phi, |_, _, _| true)?;

    // approximate the integral of fluxes along the vertical section
    let (point_ids, coords, w_along_section) = post.values_along_edges_vec(&mid_patch, &mid_section)?;
    let mut area = 0.0;
    let un = Vector::from(&[1.0, 0.0]); // it's a vertical section
    for i in 1..point_ids.len() {
        let dy = coords[i][1] - coords[i - 1][1];
        let dot_prev = vec_inner(&w_along_section[i - 1], &un);
        let dot = vec_inner(&w_along_section[i], &un);
        area += dy * (dot_prev + dot) / 2.0;
    }
    println!("Flux through vertical middle section = {}", area);
    approx_eq(area, 2.18789, 1e-5); // the analytical value is 2.1495; this is the value from Paraview

    // analytical solution
    let k = 1.0;
    let h = 12.0 - 8.0;
    let s = 0.0;
    let b = 4.0;
    let tt = 8.0;
    let s_by_tt = s / tt;
    let b_by_tt = b / tt;
    let q_by_kh = ConfinedFlow::normalized_discharge_symmetrically_placed_wall(s_by_tt, b_by_tt);
    let q_analytical = q_by_kh * k * h;
    println!("Analytical discharge: q = {}", q_analytical);

    // write Paraview files
    // let path_pvd = post.write_paraview(&mut memo, OUT_DIR, NAME)?;
    // println!("\nParaview File: {}\n", path_pvd);
    Ok(())
}
