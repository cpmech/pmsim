use crate::base::{Config, DofNumbers, Element, ParamLiquidRetention, ParamStressStrain};
use crate::StrError;
use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_tensor::Tensor2;
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct State {
    pub time: f64,
    pub delta_time: f64,
    pub primary_unknowns: Vector, // (n_equation)

    pub effective_stress: Vec<Vec<Tensor2>>, // (ncell,n_integ_point)
    pub int_values_solid: Vec<Vec<Vector>>,  // (ncell,n_integ_point)
    pub loading: Vec<Vec<bool>>,             // (ncell,n_integ_point)

    pub liquid_saturation: Vec<Vec<f64>>,    // (ncell,n_integ_point)
    pub int_values_porous: Vec<Vec<Vector>>, // (ncell,n_integ_point)
    pub wetting: Vec<Vec<bool>>,             // (ncell,n_integ_point)
}

impl State {
    pub fn new(mesh: &Mesh, dn: &DofNumbers, config: &Config) -> Result<State, StrError> {
        // gather information about element types
        let mut rods_and_beams_only = true;
        let mut has_diffusion = false;
        let mut has_solid = false;
        let mut has_porous = false;
        for cell in &mesh.cells {
            let element = dn
                .elements
                .get(&cell.attribute_id)
                .ok_or("cannot extract CellAttributeId to allocate State")?;
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

        // TODO: check compatibility of flags

        // return state with most vectors empty
        if rods_and_beams_only || has_diffusion {
            return Ok(State {
                time: 0.0,
                delta_time: 0.0,
                primary_unknowns: Vector::new(dn.n_equation),

                effective_stress: Vec::new(),
                int_values_solid: Vec::new(),
                loading: Vec::new(),

                liquid_saturation: Vec::new(),
                int_values_porous: Vec::new(),
                wetting: Vec::new(),
            });
        }

        // pre-allocate data for solid/porous elements
        let ncell = mesh.cells.len();
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
        let zero_sigma = Tensor2::new(true, mesh.ndim == 2);

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
        for cell in &mesh.cells {
            let ips = config.integ_point_data(cell)?;
            let n_integ_point = ips.len();
            let element = dn.elements.get(&cell.attribute_id).unwrap(); // already checked above
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
            time: 0.0,
            delta_time: 0.0,
            primary_unknowns: Vector::new(dn.n_equation),

            effective_stress,
            int_values_solid,
            loading,

            liquid_saturation,
            int_values_porous,
            wetting,
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::State;
    use crate::base::{Config, DofNumbers, Element, SampleParams};
    use gemlab::mesh::Samples;
    use std::collections::HashMap;

    #[test]
    fn new_works_mixed() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid_drucker_prager();
        let p3 = SampleParams::param_beam();
        let elements = HashMap::from([
            (1, Element::PorousSldLiq(p1)),
            (2, Element::Solid(p2)),
            (3, Element::Beam(p3)),
        ]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.delta_time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), ncell);
        assert_eq!(state.int_values_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.int_values_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.effective_stress[0].len(), 9 /*n_integ_point*/);
        assert_eq!(state.effective_stress[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.effective_stress[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.effective_stress[3].len(), 0 /*n_integ_point*/);
        assert_eq!(state.int_values_solid[0].len(), 0 /*n_integ_point*/);
        assert_eq!(state.int_values_solid[1].len(), 7 /*n_integ_point*/);
        assert_eq!(state.int_values_solid[2].len(), 0 /*n_integ_point*/);
        assert_eq!(state.int_values_solid[3].len(), 0 /*n_integ_point*/);
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
        let elements = HashMap::from([(1, Element::Diffusion(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), 0);
        assert_eq!(state.int_values_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.int_values_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn new_works_rod_only() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let elements = HashMap::from([(1, Element::Rod(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), 0);
        assert_eq!(state.int_values_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.int_values_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn new_works_porous_liq() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq();
        let elements = HashMap::from([(1, Element::PorousLiq(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), 0);
        assert_eq!(state.int_values_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.int_values_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_porous_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_liq_gas();
        let elements = HashMap::from([(1, Element::PorousLiqGas(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), 0);
        assert_eq!(state.int_values_solid.len(), 0);
        assert_eq!(state.loading.len(), 0);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.int_values_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = SampleParams::param_porous_sld_liq_gas();
        let elements = HashMap::from([(1, Element::PorousSldLiqGas(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), ncell);
        assert_eq!(state.int_values_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), ncell);
        assert_eq!(state.int_values_porous.len(), ncell);
        assert_eq!(state.wetting.len(), ncell);
        assert_eq!(state.liquid_saturation[0].len(), 7 /*n_integ_point*/);
    }

    #[test]
    fn new_works_solid_and_rod() {
        let mesh = Samples::mixed_shapes_2d();
        let p1 = SampleParams::param_rod();
        let p2 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Rod(p1)), (2, Element::Solid(p2))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state = State::new(&mesh, &dn, &config).unwrap();
        let ncell = mesh.cells.len();
        assert_eq!(state.time, 0.0);
        assert_eq!(state.primary_unknowns.dim(), dn.n_equation);
        assert_eq!(state.effective_stress.len(), ncell);
        assert_eq!(state.int_values_solid.len(), ncell);
        assert_eq!(state.loading.len(), ncell);
        assert_eq!(state.liquid_saturation.len(), 0);
        assert_eq!(state.int_values_porous.len(), 0);
        assert_eq!(state.wetting.len(), 0);
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_lin2();
        let p1 = SampleParams::param_rod();
        let elements = HashMap::from([(1, Element::Rod(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        let config = Config::new();
        let state_ori = State::new(&mesh, &dn, &config).unwrap();
        let state = state_ori.clone();
        let correct ="State { time: 0.0, delta_time: 0.0, primary_unknowns: NumVector { data: [0.0, 0.0, 0.0, 0.0] }, effective_stress: [], int_values_solid: [], loading: [], liquid_saturation: [], int_values_porous: [], wetting: [] }";
        assert_eq!(format!("{:?}", state), correct);
        // serialize
        let json = serde_json::to_string(&state).unwrap();
        // deserialize
        let read: State = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), correct);
    }

    #[test]
    fn new_handles_errors() {
        let mut mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let elements = HashMap::from([(1, Element::Solid(p1))]);
        let dn = DofNumbers::new(&mesh, elements).unwrap();
        mesh.cells[0].attribute_id = 100; // << never do this!
        let mut config = Config::new();
        assert_eq!(
            State::new(&mesh, &dn, &config).err(),
            Some("cannot extract CellAttributeId to allocate State")
        );
        mesh.cells[0].attribute_id = 1;
        config.n_integ_point.insert(1, 100); // wrong number
        assert_eq!(
            State::new(&mesh, &dn, &config).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }
}
