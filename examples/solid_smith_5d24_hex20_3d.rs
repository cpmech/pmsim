use gemlab::prelude::*;
use pmsim::{base::SampleMeshes, prelude::*};
use russell_lab::*;

// Smith's Example 5.24 (Figure 5.24) on page 195
//
// Smith IM, Griffiths DV, and Margetts L (2014) Programming the Finite
// Element Method, Wiley, Fifth Edition, 664p
//
// TEST GOAL
//
// This test verifies a 3D simulation with Hex20.
//
// MESH
//
// See figure on the documentation.
//
// BOUNDARY CONDITIONS
//
// Horizontally fix the vertical boundary faces perpendicular to x on the "back side" with x=0
// Horizontally fix the vertical boundary faces perpendicular to y on the "left side" with y=0
// Set all Ux,Uy,Uz to zero for the horizontal boundary faces perpendicular to z on the "bottom" with z=0
// Apply distributed load Qn = -1 on the portion of the top face with y ≤ 1
// NOTE: The "front" and "right" faces with x>0 or y>0 are NOT fixed.
//
// CONFIGURATION AND PARAMETERS
//
// Upper layer: Young = 100, Poisson = 0.3
// Lower layer: Young = 50, Poisson = 0.3
// Using reduced integration with 8 points

const NAME: &str = "ex_solid_smith_5d24_hex20_3d";
const SAVE_VTU: bool = false;

fn main() -> Result<(), StrError> {
    // mesh
    let mesh = SampleMeshes::smith_example_5d24_hex20();

    // features
    let feat = Features::new(&mesh, false);
    let faces_x_min = feat.search_faces(At::X(0.0), any_x)?;
    let faces_y_min = feat.search_faces(At::Y(0.0), any_x)?;
    let bottom = feat.search_faces(At::Z(-2.0), any_x)?;
    let top = feat.search_faces(At::Z(0.0), |x| x[1] <= 1.0)?;

    // input data
    let p1 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 100.0,
            poisson: 0.3,
        },
    };
    let p2 = ParamSolid {
        density: 1.0,
        stress_strain: ParamStressStrain::LinearElastic {
            young: 50.0,
            poisson: 0.3,
        },
    };
    let input = FemInput::new(&mesh, [(1, Element::Solid(p1)), (2, Element::Solid(p2))])?;

    // essential boundary conditions
    let zero = |_| 0.0;
    let mut essential = Essential::new();
    essential
        .on(&faces_x_min, Ebc::Ux(zero))
        .on(&faces_y_min, Ebc::Uy(zero))
        .on(&bottom, Ebc::Ux(zero))
        .on(&bottom, Ebc::Uy(zero))
        .on(&bottom, Ebc::Uz(zero));

    // natural boundary conditions
    let mut natural = Natural::new();
    natural.on(&top, Nbc::Qn(|_| -1.0));

    // configuration
    let mut config = Config::new();
    config.n_integ_point.insert(1, 8);
    config.n_integ_point.insert(2, 8);

    // FEM state
    let mut state = FemState::new(&input, &config)?;

    // FEM output
    let fn_stem = if SAVE_VTU { Some(NAME.to_string()) } else { None };
    let mut output = FemOutput::new(&input, fn_stem, None, None)?;

    // solve problem
    let mut solver = FemSolverImplicit::new(&input, &config, &essential, &natural)?;
    solver.solve(&mut state, &mut output)?;

    // check displacements
    #[rustfmt::skip]
    let uu_correct = &[
         0.000000000000000e+00,  0.000000000000000e+00, -2.246500765490787e-02,
         1.583910496060340e-03,  0.000000000000000e+00, -2.255483694735321e-02,
         3.220902825001161e-03,  0.000000000000000e+00, -2.333291718903119e-02,
         0.000000000000000e+00,  0.000000000000000e+00, -1.849079487110627e-02,
         1.543386163277156e-03,  0.000000000000000e+00, -1.883809062820367e-02,
         0.000000000000000e+00,  0.000000000000000e+00, -1.442590460790187e-02,
         7.581882957110439e-04,  0.000000000000000e+00, -1.435240682209192e-02,
         1.510921099665469e-03,  0.000000000000000e+00, -1.410677796630181e-02,
         0.000000000000000e+00,  0.000000000000000e+00, -6.163334370573628e-03,
         2.791745144538955e-03,  0.000000000000000e+00, -6.429270136326062e-03,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -2.636509763565923e-03, -2.091416880201127e-02,
         2.682119090872207e-03, -2.351562947333217e-03, -2.156639392479440e-02,
         0.000000000000000e+00,  1.920880391984773e-03, -1.258001206489784e-02,
         1.401676522438176e-03,  2.013303092837565e-03, -1.228963313721055e-02,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -4.404752518931236e-03, -1.366052652969094e-02,
         6.193924163718364e-04, -4.364483725499084e-03, -1.368326010742839e-02,
         1.280618203466245e-03, -3.937055049747129e-03, -1.386617424541054e-02,
         0.000000000000000e+00,  1.078066353252381e-03, -1.150393864686515e-02,
         8.753319133480854e-04,  1.490100951484467e-03, -1.152039549038886e-02,
         0.000000000000000e+00,  2.958188258300922e-03, -9.171834498547103e-03,
         5.428277078585079e-04,  2.965646563796275e-03, -9.093320479840279e-03,
         1.086820865799913e-03,  3.083787195165820e-03, -8.864006435062395e-03,
         0.000000000000000e+00,  2.375610447570931e-03, -4.184108581169291e-03,
         2.000018481259384e-03,  2.632220513610760e-03, -4.380144787490711e-03,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -3.579943057052639e-03, -5.979630028745818e-03,
        -7.413880481683365e-05, -3.321081990493434e-03, -5.813221407057556e-03,
         0.000000000000000e+00,  3.150650899126640e-03, -5.383111622858573e-03,
         7.612343087324263e-04,  3.312533658885393e-03, -5.104128342632701e-03,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -2.121036375171861e-03, -2.575781050442681e-03,
        -1.816710580698765e-04, -2.066534854927771e-03, -2.541288267716592e-03,
        -3.531265160052207e-04, -2.120391079723742e-03, -2.337896199085389e-03,
         0.000000000000000e+00,  2.687562242570146e-04, -2.403920299163627e-03,
         1.000591537317435e-04,  3.146936362030586e-04, -2.218104306007612e-03,
         0.000000000000000e+00,  2.651350794443137e-03, -2.133190207373987e-03,
         2.082095312685143e-04,  2.668801007310832e-03, -2.080388279569574e-03,
         4.197803511336846e-04,  2.802391949145345e-03, -1.956602157583608e-03,
         0.000000000000000e+00,  2.293393267994337e-03, -1.130507608629370e-03,
         6.752160199631226e-04,  2.510883798667245e-03, -1.197916353113878e-03,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -1.712317073186578e-03, -5.376714266281263e-04,
        -9.197039851119624e-05, -1.719430163328354e-03, -4.759289615407881e-04,
         0.000000000000000e+00,  2.044745359103121e-03, -3.707360108117185e-04,
         8.501653790242899e-05,  2.179757376749243e-03, -2.927030027983088e-04,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00, -1.707792121066834e-03,  1.459352529725203e-03,
        -3.548935758540396e-05, -1.701913033053894e-03,  1.490819285141552e-03,
        -7.238646797034549e-05, -1.669990382071488e-03,  1.451817030069681e-03,
         0.000000000000000e+00,  3.922224626300735e-04,  1.136636538839139e-03,
        -1.456332291842613e-04,  4.848030580407152e-04,  1.138391813297438e-03,
         0.000000000000000e+00,  1.668010029322225e-03,  5.260512769955115e-04,
        -9.977657940250570e-05,  1.700278490655642e-03,  5.388549291586811e-04,
        -2.111247101690338e-04,  1.825808063312428e-03,  5.433062553032318e-04,
         0.000000000000000e+00,  1.572078211332934e-03, -1.028507508147740e-04,
        -7.448485051380094e-05,  1.715786457536998e-03, -5.847263101786664e-05,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
         0.000000000000000e+00,  0.000000000000000e+00,  0.000000000000000e+00,
    ];
    vec_approx_eq(state.uu.as_data(), uu_correct, 2e-9);
    Ok(())
}
