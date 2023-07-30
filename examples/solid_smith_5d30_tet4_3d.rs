use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*, StrError};
use russell_chk::vec_approx_eq;

const FILENAME_KEY: &'static str = "ex_solid_smith_5d30_tet4_3d";

// Smith's Example 5.30 (Figure 5.30) on page 202
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a 3D simulation with Tet4.
//
// MESH
//
// See figure on the documentation.
//
// BOUNDARY CONDITIONS
//
// Horizontally fix the vertical boundary faces perpendicular to x on the "back side" with x=0
// Horizontally fix the vertical boundary faces perpendicular to y on the "left side" with y=0
// Vertically fix the horizontal boundary faces perpendicular to z on the "bottom" with z=0
// Apply vertical (Fz) concentrated loads to the top nodes:
// Fz @ 0 and 5 = -0.1667, Fz @ 1 and 4 = -0.3333
// (Do not USE more digits, as in the code, so we can compare with the Book results)
//
// CONFIGURATION AND PARAMETERS
//
// Young = 100, Poisson = 0.3

fn main() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d30_tet4();

    // features
    let find = Find::new(&mesh, None);
    let faces_x_min = find.faces(At::X(0.0), any_x)?;
    let faces_y_min = find.faces(At::Y(0.0), any_x)?;
    let bottom = find.faces(At::Z(-1.0), any_x)?;

    // parameters, DOFs, and configuration
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 100.0,
            poisson: 0.3,
        },
    };
    let data = Data::new(&mesh, [(1, Element::Solid(p1))])?;
    let config = Config::new();

    // essential boundary conditions
    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .on(&faces_x_min, Ebc::Ux(zero))
        .on(&faces_y_min, Ebc::Uy(zero))
        .on(&bottom, Ebc::Uz(zero));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural
        .at(&[0, 5], Pbc::Fz(|_| -0.1667))
        .at(&[1, 4], Pbc::Fz(|_| -0.3333));

    // simulation state
    let mut state = State::new(&data, &config)?;

    // run simulation
    let mut sim = Simulation::new(&data, &config, &essential, &natural)?;
    sim.run(&mut state)?;

    // generate Paraview file
    if false {
        let proc = PostProc::new(&mesh, &find, &data, &state);
        proc.write_vtu(&FilePath::vtu(FILENAME_KEY, true))?;
    }

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
        0.000000000000000e+00,  0.000000000000000e+00, -1.000068655710468e-02,
        2.999940775949772e-03,  0.000000000000000e+00, -9.999913345976635e-03,
        0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
        3.000109113827164e-03,  0.000000000000000e+00,  0.000000000000000e+00,
        0.000000000000000e+00,  2.999807785311942e-03, -9.999456316530526e-03,
        2.999866256528334e-03,  2.999938401332931e-03, -1.000057411788102e-02,
        0.000000000000000e+00,  3.000091440207385e-03,  0.000000000000000e+00,
        3.000108483339240e-03,  3.000132531607437e-03,  0.000000000000000e+00,
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 1e-15);
    Ok(())
}
