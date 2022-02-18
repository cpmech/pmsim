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
    use crate::{ElementConfig, SampleParams, SimConfig, StrError};
    use gemlab::mesh::Mesh;

    #[test]
    fn porous_layers_new_works() -> Result<(), StrError> {
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
        let layers = PorousLayers::new(&config)?;
        assert_eq!(layers.top_down.len(), 2);
        // assert_eq!(layers., (222, (1.0, 3.0)));
        // assert_eq!(layers[1], (111, (0.0, 1.0)));
        // assert_eq!(mesh.max[1] - mesh.min[1], layers[0].1 .1 - layers[p].1 .0);

        let mut mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let p1 = SampleParams::params_porous_sol_liq_gas(0.2, 1e-2);
        let p2 = SampleParams::params_porous_sol_liq_gas(0.4, 1e-2);
        let p3 = SampleParams::params_solid();
        config
            .elements(1, ElementConfig::PorousSolLiqGas(p1, None))?
            .elements(2, ElementConfig::PorousSolLiqGas(p2, None))?
            .elements(3, ElementConfig::Solid(p3, None))?
            .set_gravity(10.0)?; // m/s²
        let layers = PorousLayers::new(&config)?;
        assert_eq!(layers.top_down.len(), 2);
        // assert_eq!(layers[0], (2, (1.0, 3.0)));
        // assert_eq!(layers[1], (1, (0.0, 1.0)));
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
