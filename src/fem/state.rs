use super::Data;
use crate::base::{Config, Element, ParamLiquidRetention, ParamStressStrain};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

/// Holds state of a simulation, including primary and secondary variables
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct State {
    /// Time
    pub t: f64,

    /// Delta time
    pub dt: f64,

    /// Primary unknowns {U}
    ///
    /// (n_equation)
    pub uu: Vector,

    /// First time derivative of primary unknowns d{U}/dt
    ///
    /// (n_equation)
    pub vv: Vector,

    /// Second time derivative of primary unknowns d²{U}/dt²
    ///
    /// (n_equation)
    pub aa: Vector,

    /// Auxiliary time-discretization variable (Theta method)
    ///
    /// (n_equation)
    pub uu_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method)
    ///
    /// (n_equation)
    pub vv_star: Vector,

    /// Auxiliary time-discretization variable (Newmark method)
    ///
    /// (n_equation)
    pub aa_star: Vector,

    /// Total or effective stress of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub sigma: Vec<Vec<Tensor2>>,

    /// Internal values of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub ivs_solid: Vec<Vec<Vector>>,

    /// Loading flags of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub loading: Vec<Vec<bool>>,

    /// Liquid saturation of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub liquid_saturation: Vec<Vec<f64>>,

    /// Internal values for porous media of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub ivs_porous: Vec<Vec<Vector>>,

    /// wetting flags of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub wetting: Vec<Vec<bool>>,
}

impl State {
    pub fn new(data: &Data, config: &Config) -> Result<State, StrError> {
        // check number of cells
        if data.mesh.cells.len() == 0 {
            return Err("there are no cells in the mesh");
        }

        // gather information about element types
        let mut has_diffusion = false;
        let mut has_rod_or_beam = false;
        let mut has_solid = false;
        let mut has_porous_fluid = false;
        let mut has_porous_solid = false;
        for cell in &data.mesh.cells {
            let element = data.attributes.get(cell).unwrap(); // already checked by Data
            match element {
                Element::Diffusion(..) => has_diffusion = true,
                Element::Rod(..) => has_rod_or_beam = true,
                Element::Beam(..) => has_rod_or_beam = true,
                Element::Solid(..) => has_solid = true,
                Element::PorousLiq(..) => has_porous_fluid = true,
                Element::PorousLiqGas(..) => has_porous_fluid = true,
                Element::PorousSldLiq(..) => has_porous_solid = true,
                Element::PorousSldLiqGas(..) => has_porous_solid = true,
            };
        }

        // check elements
        if has_diffusion && (has_rod_or_beam || has_solid || has_porous_fluid || has_porous_solid) {
            return Err("cannot combine Diffusion elements with other elements");
        }
        if has_porous_fluid && (has_diffusion || has_rod_or_beam || has_solid || has_porous_solid) {
            return Err("cannot combine PorousLiq or PorousLiqGas with other elements");
        }

        // constants
        let ndim = data.mesh.ndim;
        let ncell = data.mesh.cells.len();
        let n_equation = data.equations.n_equation;

        // primary variables
        let t = config.control.t_ini;
        let dt = (config.control.dt)(t);
        let uu = Vector::new(n_equation);
        let (uu_star, vv, vv_star) = if config.transient || config.dynamics {
            (
                Vector::new(n_equation),
                Vector::new(n_equation),
                Vector::new(n_equation),
            )
        } else {
            (Vector::new(0), Vector::new(0), Vector::new(0))
        };
        let (aa, aa_star) = if config.dynamics {
            (Vector::new(n_equation), Vector::new(n_equation))
        } else {
            (Vector::new(0), Vector::new(0))
        };

        // return state with no secondary variables
        if !has_solid && !has_porous_fluid && !has_porous_solid {
            return Ok(State {
                t,
                dt,
                uu,
                vv,
                aa,
                uu_star,
                vv_star,
                aa_star,

                sigma: Vec::new(),
                ivs_solid: Vec::new(),
                loading: Vec::new(),

                liquid_saturation: Vec::new(),
                ivs_porous: Vec::new(),
                wetting: Vec::new(),
            });
        }

        // pre-allocate secondary arrays for solid/porous elements
        let mut sigma: Vec<Vec<Tensor2>>;
        let mut ivs_solid: Vec<Vec<Vector>>;
        let mut loading: Vec<Vec<bool>>;
        if has_solid || has_porous_solid {
            sigma = vec![Vec::new(); ncell];
            ivs_solid = vec![Vec::new(); ncell];
            loading = vec![Vec::new(); ncell];
        } else {
            sigma = Vec::new();
            ivs_solid = Vec::new();
            loading = Vec::new();
        }

        // pre-allocate secondary arrays for porous elements
        let mut liquid_saturation: Vec<Vec<f64>>;
        let ivs_porous: Vec<Vec<Vector>>;
        let mut wetting: Vec<Vec<bool>>;
        if has_porous_fluid || has_porous_solid {
            liquid_saturation = vec![Vec::new(); ncell];
            ivs_porous = vec![Vec::new(); ncell];
            wetting = vec![Vec::new(); ncell];
        } else {
            liquid_saturation = Vec::new();
            ivs_porous = Vec::new();
            wetting = Vec::new();
        }

        // function to generate integration point variables for solid/porous elements
        let two_dim = ndim == 2;
        let mut solid = |cell_id, n_integ_point, stress_strain: &ParamStressStrain| {
            // stresses
            sigma[cell_id] = (0..n_integ_point)
                .into_iter()
                .map(|_| Tensor2::new(true, two_dim))
                .collect();
            // internal values
            let n_internal_values = stress_strain.n_internal_values();
            // if n_internal_values > 0 {
            ivs_solid[cell_id] = (0..n_integ_point)
                .into_iter()
                .map(|_| Vector::new(n_internal_values))
                .collect();
            // }
            // loading flags
            // if stress_strain.elasto_plastic() {
            loading[cell_id] = vec![false; n_integ_point];
            // }
        };

        // function to generate integration point variables for porous elements
        let mut porous = |cell_id, n_integ_point, liquid_retention: &ParamLiquidRetention| {
            let sl_max = liquid_retention.max_liquid_saturation();
            liquid_saturation[cell_id] = vec![sl_max; n_integ_point];
            wetting[cell_id] = vec![false; n_integ_point];
        };

        // append (various n_integ_point) sub-vectors to secondary vectors
        for cell in &data.mesh.cells {
            let ips = config.integ_point_data(cell)?;
            let n_integ_point = ips.len();
            let element = data.attributes.get(cell).unwrap(); // already checked by Data
            match element {
                Element::Diffusion(..) => (), // unreachable because of the previous return command
                Element::Rod(..) => (),
                Element::Beam(..) => (),
                Element::Solid(p) => solid(cell.id, n_integ_point, &p.stress_strain),
                Element::PorousLiq(p) => porous(cell.id, n_integ_point, &p.retention_liquid),
                Element::PorousLiqGas(p) => porous(cell.id, n_integ_point, &p.retention_liquid),
                Element::PorousSldLiq(p) => {
                    solid(cell.id, n_integ_point, &p.stress_strain);
                    porous(cell.id, n_integ_point, &p.retention_liquid);
                }
                Element::PorousSldLiqGas(p) => {
                    solid(cell.id, n_integ_point, &p.stress_strain);
                    porous(cell.id, n_integ_point, &p.retention_liquid);
                }
            };
        }

        // initialize stresses
        if has_solid || has_porous_solid {}

        // general state
        Ok(State {
            t,
            dt,
            uu,
            vv,
            aa,
            uu_star,
            vv_star,
            aa_star,

            sigma,
            ivs_solid,
            loading,

            liquid_saturation,
            ivs_porous,
            wetting,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::State;
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::{Mesh, Samples};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_diffusion();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_rod();
        let data = Data::new(
            &mesh,
            [
                (1, Element::Diffusion(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new();
        assert_eq!(
            State::new(&data, &config).err(),
            Some("cannot combine Diffusion elements with other elements")
        );

        let p1 = SampleParams::param_porous_liq();
        let data = Data::new(
            &mesh,
            [
                (1, Element::PorousLiq(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new();
        assert_eq!(
            State::new(&data, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );

        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(
            &mesh,
            [
                (1, Element::PorousLiqGas(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Rod(p3)),
            ],
        )
        .unwrap();
        let config = Config::new();
        assert_eq!(
            State::new(&data, &config).err(),
            Some("cannot combine PorousLiq or PorousLiqGas with other elements")
        );

        let empty_mesh = Mesh {
            ndim: 2,
            points: Vec::new(),
            cells: Vec::new(),
        };
        let p1 = SampleParams::param_solid();
        let data = Data::new(&empty_mesh, [(1, Element::Solid(p1))]).unwrap();
        let config = Config::new();
        assert_eq!(State::new(&data, &config).err(), Some("there are no cells in the mesh"));

        let mesh = Samples::one_tri3();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong number
        assert_eq!(
            State::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn new_works_mixed() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid_drucker_prager();
        let p3 = SampleParams::param_beam();
        let data = Data::new(
            &mesh,
            [
                (1, Element::PorousSldLiq(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Beam(p3)),
            ],
        )
        .unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.dt, 0.1);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.sigma[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.sigma[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.sigma[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.sigma[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs_solid[0].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs_solid[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.ivs_solid[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs_solid[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[0].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.loading[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.liquid_saturation[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.liquid_saturation[1].len(), 0 /*n_integ_point*/);
        assert_eq!(state.liquid_saturation[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.liquid_saturation[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.wetting[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.wetting[1].len(), 0 /*n_integ_point*/);
        assert_eq!(state.wetting[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.wetting[3].len(), 0 /*n_integ_point*/);
    }

    #[test]
    fn new_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let mut config = Config::new();
        config.transient = true;
        let state = State::new(&data, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.vv.dim(), data.equations.n_equation);
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), data.equations.n_equation);
        assert_eq!(state.vv_star.dim(), data.equations.n_equation);
        assert_eq!(state.aa_star.dim(), 0);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.ivs_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn new_works_rod_only() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.vv.dim(), 0);
        assert_eq!(state.aa.dim(), 0);
        assert_eq!(state.uu_star.dim(), 0);
        assert_eq!(state.vv_star.dim(), 0);
        assert_eq!(state.aa_star.dim(), 0);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.ivs_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn new_works_porous_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq();
        let data = Data::new(&mesh, [(1, Element::PorousLiq(p1))]).unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_porous_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousLiqGas(p1))]).unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let data = Data::new(&mesh, [(1, Element::PorousSldLiqGas(p1))]).unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = SampleParams::param_rod();
        let p2 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Rod(p1)), (2, Element::Solid(p2))]).unwrap();
        let mut config = Config::new();
        config.dynamics = true;
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.equations.n_equation);
        assert_eq!(state.vv.dim(), data.equations.n_equation);
        assert_eq!(state.aa.dim(), data.equations.n_equation);
        assert_eq!(state.uu_star.dim(), data.equations.n_equation);
        assert_eq!(state.vv_star.dim(), data.equations.n_equation);
        assert_eq!(state.aa_star.dim(), data.equations.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.ivs_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let data = Data::new(&mesh, [(1, Element::Rod(p1))]).unwrap();
        let config = Config::new();
        let state_ori = State::new(&data, &config).unwrap();
        let state = state_ori.clone();
        let str_ori = format!("{:?}", state).to_string();
        assert!(str_ori.len() > 0);
        // serialize
        let json = serde_json::to_string(&state).unwrap();
        // deserialize
        let read: State = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), str_ori);
    }
}
