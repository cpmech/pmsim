use gemlab::prelude::*;
use plotpy::{Curve, Plot};
use pmsim::{prelude::*, StrError};

#[test]
fn test_heat_nonlinear_1d() -> Result<(), StrError> {
    // constants
    const L: f64 = 10.0;
    const T_INF: f64 = 20.0;
    const SOURCE: f64 = 5.0;
    const K_R: f64 = 2.0;
    const BETA: f64 = 0.0;

    // mesh
    const GENERATE_MESH: bool = false;
    let mesh = if GENERATE_MESH {
        let mut block = Block::new(&[[0.0, 0.0], [L, 0.0], [L, 1.0], [0.0, 1.0]])?;
        block.set_ndiv(&[10, 1])?;
        let mesh = block.subdivide(GeoKind::Qua9)?;
        mesh.write("/tmp/pmsim/mesh_heat_nonlinear_1d.dat")?;
        draw_mesh(&mesh, true, "/tmp/pmsim/mesh_heat_nonlinear_1d.svg")?;
        mesh
    } else {
        Mesh::read("data/meshes/mesh_heat_nonlinear_1d.dat")?
    };

    // features
    let find = Find::new(&mesh, None);
    let right = find.edges(At::X(L), any)?;
    let bottom = find.point_ids(At::Y(0.0), any)?;

    // parameters, DOFs, and configuration
    let p1 = ParamDiffusion {
        rho: 1.0,
        conductivity: ParamConductivity::IsotropicLinear { kr: K_R, beta: BETA },
        source: Some(SOURCE),
    };
    let data = Data::new(&mesh, [(1, Element::Diffusion(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&right, Ebc::T(|_| T_INF));

    // natural boundary conditions
    let natural = Natural::new();

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // analytical solution
    let k_inf = (1.0 + BETA * T_INF) * K_R;
    let coef = BETA * SOURCE * L * L / (2.0 * k_inf);
    let analytical = |x: f64| {
        if BETA == 0.0 {
            1.0 - f64::powf(x / L, 2.0)
        } else {
            (f64::sqrt(1.0 + 2.0 * coef * (1.0 - f64::powf(x / L, 2.0))) - 1.0) / coef
        }
    };

    // sort points by x-coordinates
    let mut pairs: Vec<_> = bottom.iter().map(|id| (*id, mesh.points[*id].coords[0])).collect();
    pairs.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());

    // extract temperatures
    let tt: Vec<_> = pairs
        .iter()
        .map(|(id, _)| state.uu[data.equations.eq(*id, Dof::T).unwrap()])
        .collect();

    // compute plot data
    let xx: Vec<_> = pairs.iter().map(|(_, x)| x / L).collect();
    let yy_num: Vec<_> = tt
        .iter()
        .map(|temp| 2.0 * k_inf * (temp - T_INF) / (SOURCE * L * L))
        .collect();
    let yy_ana: Vec<_> = pairs.iter().map(|(_, x)| analytical(*x)).collect();

    // plot
    let mut curve_num = Curve::new();
    let mut curve_ana = Curve::new();
    curve_num
        .set_line_color("#cd0000")
        .set_line_style("None")
        .set_marker_style("+");
    curve_num.draw(&xx, &yy_num);
    curve_ana.draw(&xx, &yy_ana);
    let mut plot = Plot::new();
    plot.add(&curve_ana);
    plot.add(&curve_num);
    plot.grid_and_labels("$x\\;/\\;L$", "$2\\,k_{\\infty}\\,(T - T_{\\infty})\\;/\\;(s\\,L^2)$")
        .save("/tmp/pmsim/test_heat_nonlinear_1d.svg")?;
    Ok(())
}
