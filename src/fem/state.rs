use super::Data;
use crate::base::{Config, Element, ParamLiquidRetention, ParamStressStrain};
use crate::StrError;
use russell_lab::Vector;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

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

    /// Primary unknowns {U} at the beginning of the time step
    ///
    /// (n_equation)
    pub uu_old: Vector,

    /// First time derivative d{U}/dt at the beginning of the time step
    ///
    /// (n_equation)
    pub vv_old: Vector,

    /// Second time derivative d²{U}/dt² at the beginning of the time step
    ///
    /// (n_equation)
    pub aa_old: Vector,

    /// Effective stress of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub sigma: Vec<Vec<Tensor2>>,

    /// Internal values of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub ivs: Vec<Vec<Vector>>,

    /// Loading flags of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub loading: Vec<Vec<bool>>,

    /// Liquid saturation of all cells and all integration points
    ///
    /// (ncell,n_integ_point)
    pub sl: Vec<Vec<f64>>,

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
        // gather information about element types
        let mut rods_and_beams_only = true;
        let mut has_diffusion = false;
        let mut has_solid = false;
        let mut has_porous = false;
        for cell in &data.mesh.cells {
            let element = data.attributes.get(cell).unwrap(); // already checked by Data
            match element {
                Element::Diffusion(..) => {
                    has_diffusion = true;
                }
                Element::Rod(..) => (),
                Element::Beam(..) => (),
                Element::Solid(..) => {
                    rods_and_beams_only = false;
                    has_solid = true;
                }
                Element::PorousLiq(..) => {
                    rods_and_beams_only = false;
                    has_porous = true;
                }
                Element::PorousLiqGas(..) => {
                    rods_and_beams_only = false;
                    has_porous = true;
                }
                Element::PorousSldLiq(..) => {
                    rods_and_beams_only = false;
                    has_solid = true;
                    has_porous = true;
                }
                Element::PorousSldLiqGas(..) => {
                    rods_and_beams_only = false;
                    has_solid = true;
                    has_porous = true;
                }
            };
        }

        // constants
        let ndim = data.mesh.ndim;
        let ncell = data.mesh.cells.len();
        let n_equation = data.dof_numbers.n_equation;

        // primary variables
        let t = config.control.t_ini;
        let dt = (config.control.dt)(t);
        let uu = Vector::new(n_equation);
        let (uu_old, vv, vv_old) = if config.transient || config.dynamics {
            (
                Vector::new(n_equation),
                Vector::new(n_equation),
                Vector::new(n_equation),
            )
        } else {
            (Vector::new(0), Vector::new(0), Vector::new(0))
        };
        let (aa, aa_old) = if config.dynamics {
            (Vector::new(n_equation), Vector::new(n_equation))
        } else {
            (Vector::new(0), Vector::new(0))
        };

        // TODO: check compatibility of flags

        // return state with most vectors empty
        if rods_and_beams_only || has_diffusion {
            return Ok(State {
                t,
                dt,
                uu,
                vv,
                aa,
                uu_old,
                vv_old,
                aa_old,

                sigma: Vec::new(),
                ivs: Vec::new(),
                loading: Vec::new(),

                sl: Vec::new(),
                ivs_porous: Vec::new(),
                wetting: Vec::new(),
            });
        }

        // pre-allocate data for solid/porous elements
        let mut effective_stress: Vec<Vec<Tensor2>>;
        let mut int_values_solid: Vec<Vec<Vector>>;
        let mut loading: Vec<Vec<bool>>;
        if has_solid {
            effective_stress = vec![Vec::new(); ncell];
            int_values_solid = vec![Vec::new(); ncell];
            loading = vec![Vec::new(); ncell];
        } else {
            effective_stress = Vec::new();
            int_values_solid = Vec::new();
            loading = Vec::new();
        }

        // pre-allocate data for porous elements
        let mut liquid_saturation: Vec<Vec<f64>>;
        let int_values_porous: Vec<Vec<Vector>>;
        let mut wetting: Vec<Vec<bool>>;
        if has_porous {
            liquid_saturation = vec![Vec::new(); ncell];
            int_values_porous = vec![Vec::new(); ncell];
            wetting = vec![Vec::new(); ncell];
        } else {
            liquid_saturation = Vec::new();
            int_values_porous = Vec::new();
            wetting = Vec::new();
        }

        // allocate template stress tensor
        let zero_sigma = Tensor2::new(true, ndim == 2);

        // function to generate integration point variables for solid/porous elements
        let mut solid = |cell_id, n_integ_point, stress_strain: &ParamStressStrain| {
            effective_stress[cell_id] = vec![zero_sigma.clone(); n_integ_point];
            let n_internal_values = stress_strain.n_internal_values();
            if n_internal_values > 0 {
                let zero_internal_values = Vector::new(n_internal_values);
                int_values_solid[cell_id] = vec![zero_internal_values.clone(); n_integ_point];
            }
            if stress_strain.elasto_plastic() {
                loading[cell_id] = vec![false; n_integ_point];
            }
        };

        // function to generate integration point variables for porous elements
        let mut porous = |cell_id, n_integ_point, liquid_retention: &ParamLiquidRetention| {
            let sl_max = liquid_retention.max_liquid_saturation();
            liquid_saturation[cell_id] = vec![sl_max; n_integ_point];
            wetting[cell_id] = vec![false; n_integ_point];
        };

        // update state vectors with varied sub-vectors (n_integ_point)
        for cell in &data.mesh.cells {
            let ips = config.integ_point_data(cell)?;
            let n_integ_point = ips.len();
            let element = data.attributes.get(cell).unwrap(); // already checked above
            match element {
                Element::Diffusion(..) => (), // this will not happen
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

        // general state
        Ok(State {
            t,
            dt,
            uu,
            vv,
            aa,
            uu_old,
            vv_old,
            aa_old,

            sigma: effective_stress,
            ivs: int_values_solid,
            loading,

            sl: liquid_saturation,
            ivs_porous: int_values_porous,
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
    use gemlab::mesh::Samples;

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
        assert_eq!(state.dt, 0.0);
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.sl.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.sigma[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.sigma[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.sigma[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.sigma[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs[0].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.ivs[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.ivs[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[0].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.loading[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.loading[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.sl[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.sl[1].len(), 0 /*n_integ_point*/);
        assert_eq!(state.sl[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.sl[3].len(), 0 /*n_integ_point*/);
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
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.sl.len(), 0);
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
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.sl.len(), 0);
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
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.sl.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.sl[0].len(), 7 /*n_integ_point*/);
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
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), 0);
        assert_eq!(state.ivs.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.sl.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.sl[0].len(), 7 /*n_integ_point*/);
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
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.sl.len(), ncell);
        assert_eq!(state.ivs_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.sl[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = SampleParams::param_rod();
        let p2 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Rod(p1)), (2, Element::Solid(p2))]).unwrap();
        let config = Config::new();
        let state = State::new(&data, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.t, 0.0);
        assert_eq!(state.uu.dim(), data.dof_numbers.n_equation);
        assert_eq!(state.sigma.len(), ncell);
        assert_eq!(state.ivs.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.sl.len(), 0);
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

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong number
        assert_eq!(
            State::new(&data, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }
}
