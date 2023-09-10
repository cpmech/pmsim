use gemlab::prelude::*;
use pmsim::{prelude::*, StrError};

fn main() -> Result<(), StrError> {
    // constants
    const R1: f64 = 3.0; // inner radius
    const R2: f64 = 6.0; // outer radius
    const P1: f64 = 200.0; // inner pressure (magnitude)
    const P2: f64 = 100.0; // outer pressure (magnitude)
    const YOUNG: f64 = 1000.0;
    const POISSON: f64 = 0.25;

    // filenames
    let msh_dir = "data/meshes";
    let vtu_dir = "/tmp/pmsim";
    // let msh_opt = "qua8_7701points_2500cells";
    let msh_opt = "qua8_181points_50cells";
    let fn_msh = [msh_dir, "/quarter_ring2d_", msh_opt, ".msh"].concat();
    let fn_vtu = [vtu_dir, "/quarter_ring2d_", msh_opt, ".vtu"].concat();

    // mesh
    let mesh = Mesh::from_text_file(&fn_msh)?;

    // features
    let find = Find::new(&mesh, None);
    let bottom = find.edges(At::Y(0.0), any_x)?;
    let left = find.edges(At::X(0.0), any_x)?;
    let inner_circle = find.edges(At::Circle(0.0, 0.0, 3.0), any_x)?;
    let outer_circle = find.edges(At::Circle(0.0, 0.0, 6.0), any_x)?;

    for f in &inner_circle {
        println!("{:?}", f);
    }
    // println!("{:?}", outer_circle);

    // parameters, DOFs, and configuration
    let param1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: YOUNG,
            poisson: POISSON,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(param1))])?;
    let config = Config::new();

    // essential boundary conditions
    let mut essential = Essential::new();
    essential.on(&left, Ebc::Ux(|_| 0.0)).on(&bottom, Ebc::Uy(|_| 0.0));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .on(&inner_circle, Nbc::Qn(|_| -P1))
        .on(&outer_circle, Nbc::Qn(|_| -P2));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // analytical solution: radial displacement (ur)
    // Reference (page 160)
    // Sadd MH (2005) Elasticity: Theory, Applications and Numerics, Elsevier, 474p
    let rr1 = R1 * R1;
    let rr2 = R2 * R2;
    let drr = rr2 - rr1;
    let dp = P2 - P1;
    let aa = rr1 * rr2 * dp / drr;
    let bb = (rr1 * P1 - rr2 * P2) / drr;
    let c1 = (1.0 + POISSON) / YOUNG;
    let c2 = 1.0 - 2.0 * POISSON;
    let analytical_ur = |r: f64| c1 * (r * c2 * bb - aa / r);

    // results
    println!("ur(3.0) = {:?}", analytical_ur(3.0));
    println!("ur(6.0) = {:?}", analytical_ur(6.0));

    // check displacements at bottom
    println!("");
    let mut selection = find.point_ids(At::Y(0.0), any_x)?;
    selection.sort_by(|a, b| {
        mesh.points[*a].coords[0]
            .partial_cmp(&mesh.points[*b].coords[0])
            .unwrap()
    });
    for p in &selection {
        let r = mesh.points[*p].coords[0];
        let eq = data.equations.eq(*p, Dof::Ux).unwrap();
        let ux = state.uu[eq];
        let diff = f64::abs(ux - analytical_ur(r));
        println!("{:3}, r = {:.3}, Ux = {:18.15}, diff = {:.2e}", p, r, ux, diff);
        assert!(diff < 1e-5);
    }

    // check displacements at left
    println!("");
    let mut selection = find.point_ids(At::X(0.0), any_x)?;
    selection.sort_by(|a, b| {
        mesh.points[*a].coords[1]
            .partial_cmp(&mesh.points[*b].coords[1])
            .unwrap()
    });
    for p in &selection {
        let r = mesh.points[*p].coords[1];
        let eq = data.equations.eq(*p, Dof::Uy).unwrap();
        let uy = state.uu[eq];
        let diff = f64::abs(uy - analytical_ur(r));
        println!("{:3}, r = {:.3}, Ux = {:18.15}, diff = {:.2e}", p, r, uy, diff);
        assert!(diff < 1e-5);
    }

    // post-processing
    let post = PostProc::new(&mesh, &find, &data, &state);
    post.write_vtu(&fn_vtu)?;
    Ok(())
}
