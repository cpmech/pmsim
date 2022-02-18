#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ElementConfig, ModelPorousSolLiq, ModelPorousSolLiqGas, ModelRealDensity, SimConfig, StrError};
use gemlab::mesh::{At, CellAttributeId, CellId, PointId};
use russell_tensor::Tensor2;
use std::collections::{HashMap, HashSet};

/// Holds (y_min,y_max) for 2D or (z_min,z_max) for 3D
type MinMaxElevation = (f64, f64);

/// Holds the (attribute,limits) pair for a layer
type LayerLimits = (CellAttributeId, MinMaxElevation);

/// Holds essential information about a porous layer which corresponds to a CellAttributeId
#[derive(Clone, Copy, Debug)]
pub struct PorousLayer {
    pub id: CellAttributeId, // identification number
    pub z_min: f64,          // minimum elevation (y in 2D or z in 3D)
    pub z_max: f64,          // maximum elevation (y in 2D or z in 3D)
    pub overburden: f64,     // vertical stress (total, **not** effective) at the top (z_max) of the layer
}

/// Holds layers in a supposedly rectangular (2D) or parallepiped (3D) mesh
///
/// # Note
///
/// * The layers are found by looking at the cells near the origin,
///   thus you have to make sure that the mesh has horizontal layers
///   with all cells in the same layer having the same CellAttributeId
#[derive(Clone, Debug)]
pub struct PorousLayers {
    pub top_down: Vec<PorousLayer>, // layers sorted from top to down
}

impl PorousLayers {
    /// Finds layers by looking at the cells near the origin
    ///
    /// # Note
    ///
    /// * The overburden stress is set to zero and must be computed later on
    pub fn new(config: &SimConfig) -> Result<Self, StrError> {
        // mesh and space_ndim
        let mesh = config.mesh;
        let space_ndim = mesh.space_ndim;

        // find cells near x_min
        let mut cells_near_x_min: HashSet<CellId> = HashSet::new();
        if mesh.space_ndim == 2 {
            let edge_keys = mesh.find_boundary_edges(At::X(mesh.coords_min[0]))?;
            if edge_keys.len() < 1 {
                return Err("cannot find at least one vertical edge at x_min");
            }
            for edge_key in &edge_keys {
                let edge = mesh.boundary_edges.get(edge_key).unwrap();
                for cell_id in &edge.shared_by_2d_cells {
                    cells_near_x_min.insert(*cell_id);
                }
            }
        } else {
            let face_keys = mesh.find_boundary_faces(At::X(mesh.coords_min[0]))?;
            if face_keys.len() < 1 {
                return Err("cannot find at least one vertical face at x_min");
            }
            for face_key in &face_keys {
                let face = mesh.boundary_faces.get(face_key).unwrap();
                for cell_id in &face.shared_by_cells {
                    cells_near_x_min.insert(*cell_id);
                }
            }
        }

        // find layers (unsorted)
        let mut layers: HashMap<CellAttributeId, PorousLayer> = HashMap::new();
        for cell_id in &cells_near_x_min {
            let cell = &mesh.cells[*cell_id];
            let attribute_id = cell.attribute_id;
            let element_config = config
                .element_configs
                .get(&attribute_id)
                .ok_or("cannot find CellAttributeId in SimConfig")?;
            match element_config {
                ElementConfig::PorousSolLiq(..) | ElementConfig::PorousSolLiqGas(..) => (),
                _ => continue, // skip other elements
            }
            let shape = mesh.alloc_shape_cell(*cell_id)?;
            match layers.get_mut(&attribute_id) {
                Some(layer) => {
                    layer.z_min = f64::min(layer.z_min, shape.coords_min[space_ndim - 1]);
                    layer.z_max = f64::max(layer.z_max, shape.coords_max[space_ndim - 1]);
                }
                None => {
                    layers.insert(
                        attribute_id,
                        PorousLayer {
                            id: attribute_id,
                            z_min: shape.coords_min[space_ndim - 1],
                            z_max: shape.coords_max[space_ndim - 1],
                            overburden: 0.0,
                        },
                    );
                }
            };
        }

        // convert map to vec and sort
        let mut top_down: Vec<_> = layers.keys().map(|k| *layers.get(k).unwrap()).collect();
        top_down.sort_by(|a, b| b.z_min.partial_cmp(&a.z_min).unwrap());
        if top_down.len() < 1 {
            return Err("at least one layer (CellAttributeId) is required");
        }

        // total height of the full set of layers
        let height = top_down.first().unwrap().z_max - top_down.last().unwrap().z_min;
        assert!(height > 0.0);

        // compute overburden stress
        let gravity = config.gravity;
        let two_dim = space_ndim == 2;
        let mut cumulated_overburden_stress = 0.0;
        for layer in &mut top_down {
            layer.overburden = cumulated_overburden_stress; // the first layer at the top has zero overburden
            let thickness = layer.z_max - layer.z_min;
            assert!(thickness > 0.0);
            let base_elevation = layer.z_min;
            let element_config = config.element_configs.get(&layer.id).unwrap();
            let rho_ini = match element_config {
                ElementConfig::PorousSolLiq(params, _) => {
                    let porous = ModelPorousSolLiq::new(params, two_dim)?;
                    let pl = porous
                        .model_density_liquid
                        .pressure_at_elevation(base_elevation, height, gravity)?;
                    porous.calc_rho_ini(pl)?
                }
                ElementConfig::PorousSolLiqGas(params, _) => {
                    let porous = ModelPorousSolLiqGas::new(params, two_dim)?;
                    let pl = porous
                        .model_density_liquid
                        .pressure_at_elevation(base_elevation, height, gravity)?;
                    let pg = porous
                        .model_density_gas
                        .pressure_at_elevation(base_elevation, height, gravity)?;
                    porous.calc_rho_ini(pl, pg)?
                }
                _ => panic!("INTERNAL ERROR: only porous models should be in the layers vector"),
            };
            let delta_sigma_v = rho_ini * config.gravity * thickness;
            cumulated_overburden_stress += delta_sigma_v;
        }

        // done
        Ok(PorousLayers { top_down })
    }
}

/// Implements geostatic stress state calculator
#[derive(Debug)]
pub struct Geostatics {
    layers: PorousLayers,
}

impl Geostatics {
    /// Returns a new StateGeostatic instance
    ///
    /// # Note
    ///
    /// * Geostatics initialization requires a rectangular (2D) or parallepiped (3D) mesh
    /// * The edges or faces at the minimum coordinates forming a vertical column are used to determine the layers
    /// * The datum is at y=y_min (2D) or z=z_min (3D)
    /// * The water table is at y=y_max (2D) or z=z_max (3D), thus only fully water-saturated (with sl_max) states are considered
    pub fn new(config: &SimConfig) -> Result<Self, StrError> {
        let mesh = config.mesh;
        let space_ndim = mesh.space_ndim;
        if space_ndim < 2 || space_ndim > 3 {
            return Err("geostatics initialization requires space_ndim = 2 or 3");
        }
        Ok(Geostatics {
            layers: PorousLayers::new(config)?,
        })
    }

    // pub fn calc_liquid_pressure(&self, coords: &[f64]) -> Result<f64, StrError> {
    //     Ok(0.0)
    // }

    // Calculates effective stresses, liquid pressure, and gas pressure
    // pub fn calc_stress(&self, _elevation: f64) -> Result<(Tensor2, f64, f64), StrError> {
    //     let stress_effective = Tensor2::new(true, self.config.two_dim);
    //     let (p_l, p_g) = (0.0, 0.0);
    //     Ok((stress_effective, p_l, p_g))
    // }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Geostatics, PorousLayers};
    use crate::{ElementConfig, ParamPorousSolLiqGas, ParamSolid, SampleParams, SimConfig, StrError};
    use gemlab::mesh::Mesh;
    use russell_chk::assert_approx_eq;

    // Returns the parameters of a two-layers column for testing
    fn two_layers_slg(height: f64, z_middle: f64) -> (ParamSolid, ParamPorousSolLiqGas, ParamPorousSolLiqGas, f64) {
        // parameters
        let footing = SampleParams::params_solid();
        let upper = SampleParams::params_porous_sol_liq_gas(0.4, 1e-2);
        let lower = SampleParams::params_porous_sol_liq_gas(0.2, 1e-2);
        // constants
        let g = 10.0; // gravity
        let sl_max = 0.95; // max liquid saturation
        let hh = height; // height of porous column
        let z = z_middle; // at bottom of upper layer
        let nf0 = upper.porosity_initial;
        let rho_ss = upper.density_solid;
        // liquid
        let rho_l_ref = upper.density_liquid.rho_ref;
        let pl_ref = upper.density_liquid.p_ref;
        let cc_l = upper.density_liquid.cc;
        let pl_mid = pl_ref + (rho_l_ref / cc_l) * (f64::exp(g * cc_l * (hh - z)) - 1.0);
        let rho_l_mid = rho_l_ref + cc_l * (pl_mid - pl_ref);
        // println!("pl_bottom = {}, rho_l_bottom = {}", pl_bottom, rho_l_bottom);
        assert_approx_eq!(pl_mid, 2.0 * g * 1.0, 1e-4);
        assert_approx_eq!(rho_l_mid, 1.0, 1e-5);
        // gas
        let rho_g_ref = upper.density_gas.rho_ref;
        let pg_ref = upper.density_gas.p_ref;
        let cc_g = upper.density_gas.cc;
        let pg_bottom = pg_ref + (rho_g_ref / cc_g) * (f64::exp(g * cc_g * (hh - z)) - 1.0);
        let rho_g_bottom = rho_g_ref + cc_g * (pg_bottom - pg_ref);
        // println!("pg_bottom = {}, rho_g_bottom = {}", pg_bottom, rho_g_bottom);
        assert_approx_eq!(pg_bottom, 2.0 * g * 0.0012, 1e-5);
        assert_approx_eq!(rho_g_bottom, 0.0012, 1e-6);
        // mixture
        let rho_mid = (1.0 - nf0) * rho_ss + nf0 * sl_max * rho_l_mid + nf0 * (1.0 - sl_max) * rho_g_bottom;
        // println!("rho_bottom = {}", rho_bottom);
        // done
        (footing, upper, lower, rho_mid)
    }

    #[test]
    fn porous_layers_new_works() -> Result<(), StrError> {
        let (footing, upper, lower, rho_mid) = two_layers_slg(3.0, 1.0);
        let mut mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        config
            .elements(111, ElementConfig::PorousSolLiqGas(lower, None))?
            .elements(222, ElementConfig::PorousSolLiqGas(upper, None))?
            .elements(333, ElementConfig::Solid(footing, None))?
            .set_gravity(10.0)?; // m/s²

        let layers = PorousLayers::new(&config)?;
        assert_eq!(layers.top_down.len(), 2);
        let top = &layers.top_down[0];
        assert_eq!(top.id, 222);
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.overburden, 0.0);
        let bottom = &layers.top_down[1];
        assert_eq!(bottom.id, 111);
        assert_eq!(bottom.z_min, 0.0);
        assert_eq!(bottom.z_max, 1.0);
        assert_approx_eq!(bottom.overburden, 2.0 * 10.0 * rho_mid, 1e-15);

        let mut mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let p1 = SampleParams::params_porous_sol_liq_gas(0.2, 1e-2);
        let p2 = SampleParams::params_porous_sol_liq_gas(0.4, 1e-2);
        let p3 = SampleParams::params_solid();
        config
            .elements(1, ElementConfig::PorousSolLiqGas(lower, None))?
            .elements(2, ElementConfig::PorousSolLiqGas(upper, None))?
            .elements(3, ElementConfig::Solid(footing, None))?
            .set_gravity(10.0)?; // m/s²
        let layers = PorousLayers::new(&config)?;
        assert_eq!(layers.top_down.len(), 2);
        let top = &layers.top_down[0];
        assert_eq!(top.id, 2);
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.overburden, 0.0);
        let bottom = &layers.top_down[1];
        assert_eq!(bottom.id, 1);
        assert_eq!(bottom.z_min, 0.0);
        assert_eq!(bottom.z_max, 1.0);
        assert_approx_eq!(bottom.overburden, 2.0 * 10.0 * rho_mid, 1e-15);
        Ok(())
    }

    #[test]
    fn geostatics_new_works() -> Result<(), StrError> {
        let mut mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let p111 = SampleParams::params_porous_sol_liq_gas(0.2, 1e-2);
        let p222 = SampleParams::params_porous_sol_liq_gas(0.4, 1e-2);
        let p333 = SampleParams::params_solid();
        config
            .elements(111, ElementConfig::PorousSolLiqGas(p111, None))?
            .elements(222, ElementConfig::PorousSolLiqGas(p222, None))?
            .elements(333, ElementConfig::Solid(p333, None))?
            .set_gravity(10.0)?; // m/s²
        Geostatics::new(&config)?;
        Ok(())
    }
}
