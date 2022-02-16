#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ElementConfig, ModelRealDensity, SimConfig, StrError};
use gemlab::mesh::{At, CellAttributeId, CellId, PointId};
use russell_tensor::Tensor2;
use std::collections::{HashMap, HashSet};

/// Holds information about a porous layer to be used in geostatics computations
///
/// # Note
///
/// This struct assists in calculating densities and pressures along a porous column.
/// The code considers **maximum liquid saturation** and minimum gas saturation.
///
/// * `rho_l_real` -- is the intrinsic (real) liquid density `ρL`
/// * `rho_g_real` -- is the intrinsic (real) gas density `ρG`
/// * `sat_liq` -- is the liquid saturation, e.g., 1.0 or 0.95 `sl`
/// * `sat_gas` -- is the gas saturation `sg = 1 - sl`
/// * `rho_mixture` -- is the partial density of the porous medium (mixture) `ρ = nf・sl・ρL + nf・sg・ρG + (1-nf)・ρS`
/// * `n_s` -- is the partial fraction of solids `ns`
/// * `sig_v_total` -- is the total vertical stress at an elevation `z` with `H` being the ... `σV_total  = σV_0 + ρ・g・(H - z)`
///
struct PorousLayer {
    // geometry
    point_ids: Vec<PointId>,             // ids of points in this layer
    attribute_ids: Vec<CellAttributeId>, // attribute ids of cells within this layer
    z_min: f64,                          // coordinate (elevation) at the bottom of this layer
    z_max: f64,                          // coordinate (elevation) at the top of this layer

    // parameters: porous medium
    sl_max: f64,       // maximum liquid saturation; e.g. 1.0
    nf_0: f64,         // (nf0) initial (constant) porosity
    rho_s_real_0: f64, // (rhoS0) initial density of solids

    // parameters: total stress analysis
    tot_rho: f64,     // density for total stress analyses
    tot_stress: bool, // total stress analysis

    // additional data
    kk0: f64, // coefficient to multiply effective vertical stresses and obtain horizontal effective stresses
    sig_v_total_top: f64, // state @ top of layer

    // auxiliary
    height: f64,                 // maximum height
    gravity: f64,                // gravity
    model_liq: ModelRealDensity, // liquid model
    model_gas: ModelRealDensity, // gas model
}

impl PorousLayer {
    // Calculates the liquid real density (rho_l_real) and pressure (pl) at given elevation
    pub fn calc_liquid_density_and_pressure(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_l_real = self
            .model_liq
            .density_at_elevation(elevation, self.height, self.gravity)?;
        let pl = self
            .model_liq
            .pressure_at_elevation(elevation, self.height, self.gravity)?;
        Ok((rho_l_real, pl))
    }

    // Calculates the gas real density (rho_g_real) and pressure (pg) at given elevation
    pub fn calc_gas_density_and_pressure(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_g_real = self
            .model_gas
            .density_at_elevation(elevation, self.height, self.gravity)?;
        let pg = self
            .model_gas
            .pressure_at_elevation(elevation, self.height, self.gravity)?;
        Ok((rho_g_real, pg))
    }

    // Calculates the porous medium partial density (rho_mixture) and total vertical stress (sigma_v_total)
    pub fn calc_mixture_density_and_total_vertical_stress(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        let rho_mixture = if self.tot_stress {
            self.tot_rho
        } else {
            let (rho_l_real, pl) = self.calc_liquid_density_and_pressure(elevation)?;
            let (rho_g_real, pg) = self.calc_gas_density_and_pressure(elevation)?;
            let sl = self.sl_max;
            let sg = 1.0 - sl;
            let nf = self.nf_0;
            let rho_s_real = self.rho_s_real_0;
            nf * sl * rho_l_real + nf * sg * rho_g_real + (1.0 - nf) * rho_s_real
        };
        let delta_z = self.z_max - elevation;
        let sigma_v_total = self.sig_v_total_top + rho_mixture * self.gravity * delta_z;
        Ok((rho_mixture, sigma_v_total))
    }
}

type MinMaxElevation = (f64, f64);

type LayerLimits = (CellAttributeId, MinMaxElevation);

fn find_layers_top_to_bottom(config: &SimConfig) -> Result<Vec<LayerLimits>, StrError> {
    // mesh and space_ndim
    let mesh = config.mesh;
    let space_ndim = mesh.space_ndim;

    // find cells near x_min
    let mut cells_near_x_min: HashSet<CellId> = HashSet::new();
    if mesh.space_ndim == 2 {
        let edge_keys = mesh.find_boundary_edges(At::X(mesh.min[0]))?;
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
        let face_keys = mesh.find_boundary_faces(At::X(mesh.min[0]))?;
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
    let mut layers: HashMap<CellAttributeId, MinMaxElevation> = HashMap::new();
    for cell_id in &cells_near_x_min {
        let cell = &mesh.cells[*cell_id];
        let attribute_id = cell.attribute_id;
        let element_config = config
            .element_configs
            .get(&attribute_id)
            .ok_or("cannot find CellAttributeId in SimConfig")?;
        match &element_config {
            ElementConfig::PorousSolLiq(..) | ElementConfig::PorousSolLiqGas(..) => (),
            _ => continue, // skip other elements
        }
        let shape = mesh.alloc_shape_cell(*cell_id)?;
        match layers.get_mut(&attribute_id) {
            Some(limits) => {
                limits.0 = f64::min(limits.0, shape.min_coords[space_ndim - 1]);
                limits.1 = f64::max(limits.1, shape.max_coords[space_ndim - 1]);
            }
            None => {
                layers.insert(
                    attribute_id,
                    (shape.min_coords[space_ndim - 1], shape.max_coords[space_ndim - 1]),
                );
            }
        };
    }

    // convert map to vec
    let mut layers_top_to_bottom: Vec<_> = layers.keys().map(|k| (*k, *layers.get(k).unwrap())).collect();

    // sort layers top-to-bottom
    layers_top_to_bottom.sort_by(|(_, a_limits), (_, b_limits)| b_limits.0.partial_cmp(&a_limits.0).unwrap());
    Ok(layers_top_to_bottom)
}

/// Implements geostatic stress state calculator
#[derive(Debug)]
pub struct Geostatics {}

impl Geostatics {
    /// Returns a new StateGeostatic instance
    ///
    /// # Note
    ///
    /// * Geostatics initialization requires a rectangular (parallepiped) mesh
    /// * The edges or faces at the minimum coordinates forming a vertical column are used to determine the layers
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max (2D) or z=z_max (3D), thus only fully water-saturated states are considered
    pub fn new(config: &SimConfig) -> Result<Self, StrError> {
        // check
        let mesh = config.mesh;
        let space_ndim = mesh.space_ndim;
        if space_ndim < 2 || space_ndim > 3 {
            return Err("geostatics initialization requires space_ndim = 2 or 3");
        }

        // find layers
        let layers_top_to_bottom = find_layers_top_to_bottom(config)?;

        // compute overburden stress
        let mut cumulated_overburden = 0.0;
        for (attribute_id, limits) in &layers_top_to_bottom {
            let element_config = config.element_configs.get(attribute_id).unwrap();
            println!("{:?}", element_config);
            // TODO
        }

        Ok(Geostatics {})
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
    use super::{find_layers_top_to_bottom, Geostatics};
    use crate::{ElementConfig, SampleParams, SimConfig, StrError};
    use gemlab::mesh::Mesh;

    #[test]
    fn find_layers_top_to_bottom_works() -> Result<(), StrError> {
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
        let layers = find_layers_top_to_bottom(&config)?;
        assert_eq!(layers[0], (222, (1.0, 3.0)));
        assert_eq!(layers[1], (111, (0.0, 1.0)));

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
        let layers = find_layers_top_to_bottom(&config)?;
        assert_eq!(layers[0], (2, (1.0, 3.0)));
        assert_eq!(layers[1], (1, (0.0, 1.0)));
        Ok(())
    }

    #[test]
    fn new_works() -> Result<(), StrError> {
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
