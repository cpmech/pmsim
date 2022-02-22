use crate::{ElementConfig, ModelPorous, SimConfig, StrError};
use gemlab::mesh::{At, CellAttributeId, CellId};
use russell_tensor::Tensor2;
use std::collections::{HashMap, HashSet};

/// Holds essential information about a porous layer
#[derive(Clone, Copy)]
struct LayerInfo {
    attribute_id: CellAttributeId, // identification number = CellAttributeId
    z_min: f64,                    // minimum elevation (y in 2D or z in 3D)
    z_max: f64,                    // maximum elevation (y in 2D or z in 3D)
}

/// Holds data for a porous layer corresponding to a CellAttributeId
pub struct PorousLayer {
    /// identification number = CellAttributeId
    pub attribute_id: CellAttributeId,

    /// minimum elevation (y in 2D or z in 3D)
    pub z_min: f64,

    /// maximum elevation (y in 2D or z in 3D)
    pub z_max: f64,

    /// total (**not** effective) vertical stress at the top (z_max) of the layer
    /// positive values means compression (**soil** mechanics convention)
    pub overburden: f64,

    /// material models
    pub model: ModelPorous,
}

/// Implements geostatics for a rectangular (2D) or parallepiped (3D) mesh
///
/// # Note
///
/// * This requires a rectangular (2D) or parallepiped (3D) mesh (not verified)
/// * The layers are found by looking at the cells near the origin,
///   thus you have to make sure that the mesh has horizontal layers
///   with all cells in the same layer having the same CellAttributeId
/// * The edges or faces at the minimum coordinates forming a vertical column
///   are used to determine the layers
/// * The datum is at y=y_min (2D) or z=z_min (3D)
/// * The water table is at y=y_max (2D) or z=z_max (3D),
///   thus only **fully liquid-saturated** (with sl_max) states are considered
pub struct Geostatics {
    /// number of dimensions
    pub space_ndim: usize,

    /// gravity constant
    pub gravity: f64,

    /// the height of the porous domain (column) = height of the whole set of layers
    pub height: f64,

    /// layers sorted from top to bottom
    pub layers: Vec<PorousLayer>,
}

impl Geostatics {
    /// Allocate a new instance
    ///
    /// Also discovers the layers from the Mesh, computes layers
    /// limits and allocates models.
    pub fn new(config: &SimConfig) -> Result<Self, StrError> {
        // mesh and space_ndim
        let mesh = config.mesh;
        let space_ndim = mesh.space_ndim;
        let two_dim = space_ndim == 2;
        if space_ndim < 2 || space_ndim > 3 {
            return Err("Geostatics requires space_ndim = 2 or 3");
        }

        // param for fluids
        let param_fluids = match &config.param_fluids {
            Some(p) => p,
            None => return Err("param for fluids must be set first"),
        };

        // find cells near x_min
        let x_min = mesh.coords_min[0];
        let mut cells_near_x_min: HashSet<CellId> = HashSet::new();
        if space_ndim == 2 {
            let edge_keys = mesh.find_boundary_edges(At::X(x_min))?;
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
            let face_keys = mesh.find_boundary_faces(At::X(x_min))?;
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

        // collect attribute ids and z_min/z_max of layers
        let mut infos: HashMap<CellAttributeId, LayerInfo> = HashMap::new();
        for cell_id in &cells_near_x_min {
            let cell = &mesh.cells[*cell_id];
            let element_config = config.get_element_config(cell.attribute_id)?;
            match element_config {
                ElementConfig::Porous(..) => (),
                _ => continue, // skip other types
            };
            let shape = mesh.alloc_shape_cell(*cell_id)?;
            match infos.get_mut(&cell.attribute_id) {
                Some(info) => {
                    info.z_min = f64::min(info.z_min, shape.coords_min[space_ndim - 1]);
                    info.z_max = f64::max(info.z_max, shape.coords_max[space_ndim - 1]);
                }
                None => {
                    infos.insert(
                        cell.attribute_id,
                        LayerInfo {
                            attribute_id: cell.attribute_id,
                            z_min: shape.coords_min[space_ndim - 1],
                            z_max: shape.coords_max[space_ndim - 1],
                        },
                    );
                }
            };
        }

        // convert map to vec, then sort
        let mut top_down: Vec<_> = infos.keys().map(|k| *infos.get(k).unwrap()).collect();
        top_down.sort_by(|a, b| b.z_min.partial_cmp(&a.z_min).unwrap());
        if top_down.len() < 1 {
            return Err("at least one layer (CellAttributeId) is required");
        }

        // allocate vector of porous layers
        let mut layers: Vec<PorousLayer> = Vec::new();
        for info in &top_down {
            let element_config = config.get_element_config(info.attribute_id)?;
            let model = match element_config {
                ElementConfig::Porous(param_porous, _) => ModelPorous::new(param_fluids, param_porous, two_dim)?,
                _ => panic!("INTERNAL ERROR: element_config is missing"), // not supposed to happen
            };
            layers.push(PorousLayer {
                attribute_id: info.attribute_id,
                z_min: info.z_min,
                z_max: info.z_max,
                overburden: 0.0,
                model,
            })
        }

        // total height of the whole set of layers
        let height = layers.first().unwrap().z_max - layers.last().unwrap().z_min;
        assert!(height > 0.0);

        // compute overburden stress
        let gravity = config.gravity;
        let mut cumulated_overburden_stress = 0.0;
        for layer in &mut layers {
            layer.overburden = cumulated_overburden_stress; // the first layer at the top has zero overburden
            let delta_sigma_v =
                layer
                    .model
                    .calc_sigma_z_ini(layer.overburden, layer.z_min, height, layer.z_max, gravity)?;
            cumulated_overburden_stress += delta_sigma_v;
        }

        // done
        Ok(Geostatics {
            space_ndim,
            gravity,
            height,
            layers,
        })
    }

    /// Calculates liquid pressure at given elevation
    pub fn calc_pl(&self, elevation: f64) -> Result<f64, StrError> {
        for layer in &self.layers {
            if elevation >= layer.z_min && elevation <= layer.z_max {
                return layer.model.calc_pl(elevation, self.height, self.gravity);
            }
        }
        Err("elevation is outside the porous region limits")
    }

    /// Calculates liquid and gas pressure at given elevation
    pub fn calc_pl_and_pg(&self, elevation: f64) -> Result<(f64, f64), StrError> {
        for layer in &self.layers {
            if elevation >= layer.z_min && elevation <= layer.z_max {
                let pl = layer.model.calc_pl(elevation, self.height, self.gravity)?;
                let pg = layer.model.calc_pg(elevation, self.height, self.gravity)?;
                return Ok((pl, pg));
            }
        }
        Err("elevation is outside the porous region limits")
    }

    /// Calculates effective stress at given elevation (e.g., integration point)
    ///
    /// Note: negative stress component means compression according to continuum/solid mechanics
    pub fn calc_effective_stress(&self, elevation: f64) -> Result<Tensor2, StrError> {
        let (symmetric, two_dim) = (true, self.space_ndim == 2);
        for layer in &self.layers {
            if elevation >= layer.z_min && elevation <= layer.z_max {
                let pl = layer.model.calc_pl(elevation, self.height, self.gravity)?;
                let rho_ini = 0.0; // TODO layer.model.calc_rho_ini(elevation, self.height, self.gravity)?;
                let sigma_v_total = layer.overburden + rho_ini * self.gravity * (layer.z_max - elevation);
                let sigma_v_effective = sigma_v_total - pl;
                let sigma_h_effective = layer.model.kk0 * sigma_v_effective;
                let (sig_x, sig_y, sig_z) = if two_dim {
                    (-sigma_h_effective, -sigma_v_effective, -sigma_h_effective)
                } else {
                    (-sigma_h_effective, -sigma_h_effective, -sigma_v_effective)
                };
                let mut stress = Tensor2::new(symmetric, two_dim);
                stress.vec[0] = sig_x;
                stress.vec[1] = sig_y;
                stress.vec[2] = sig_z;
                return Ok(stress);
            }
        }
        Err("elevation is outside the porous region limits")
    }

    /// Finds the cell attribute id of the layer containing a given elevation
    pub fn find_attribute_id(&self, elevation: f64) -> Result<CellAttributeId, StrError> {
        for layer in &self.layers {
            if elevation >= layer.z_min && elevation <= layer.z_max {
                return Ok(layer.attribute_id);
            }
        }
        Err("elevation is outside the porous region limits")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Geostatics;
    use crate::{ElementConfig, ParamPorous, ParamSolid, SampleParam, SimConfig, StrError};
    use gemlab::mesh::Mesh;
    use russell_chk::assert_approx_eq;

    // Returns the parameters of a two-layer porous column with a footing at the top and
    // the total vertical stress at the "middle" (transition between upper and lower layers)
    // Returns: (footing, upper, lower, sigma_v_mid)
    fn two_layers(
        height: f64,
        z_middle: f64,
        with_gas: bool,
        incompressible_liq: bool,
    ) -> (ParamSolid, ParamPorous, ParamPorous, f64) {
        // constants
        let hh = height; // height of porous column
        let z = z_middle; // at the middle of upper layer
        let g = 10.0; // gravity
        let nf0 = 0.4; // initial porosity
        let rho_ss = 2.7; // initial/constant real density of solids

        // auxiliary constants
        let gg_s = rho_ss; // specific gravity of solids
        let e0 = nf0 / (1.0 - nf0); // initial void ratio

        // constants: liquid
        let sl_max = if with_gas { 0.95 } else { 1.0 }; // max liquid saturation
        let rho_l_ref = 1.0; // reference liquid density
        let pl_ref = 0.0; // reference liquid pressure
        let cc_l = if incompressible_liq { 1e-12 } else { 4.53e-7 }; // liquid compressibility

        // constants: gas
        let rho_g_ref = if with_gas { 0.0012 } else { 0.0 }; // reference gas density
        let pg_ref = 0.0; // reference gas pressure
        let cc_g = 1.17e-5; // gas compressibility

        // pressure and density: liquid
        let pl_mid = pl_ref + (rho_l_ref / cc_l) * (f64::exp(g * cc_l * (hh - z)) - 1.0);
        let rho_l_mid = rho_l_ref + cc_l * (pl_mid - pl_ref);
        assert_approx_eq!(pl_mid, (hh - z) * g * rho_l_ref, 1e-4);
        assert_approx_eq!(rho_l_mid, rho_l_ref, 1e-5);

        // pressure and density: gas
        let pg_mid = pg_ref + (rho_g_ref / cc_g) * (f64::exp(g * cc_g * (hh - z)) - 1.0);
        let rho_g_mid = rho_g_ref + cc_g * (pg_mid - pg_ref);
        assert_approx_eq!(pg_mid, (hh - z) * g * rho_g_ref, 1e-5);
        assert_approx_eq!(rho_g_mid, rho_g_ref, 1e-6);

        // partial density of the mixture and total vertical stress at the middle
        let rho_mid = (1.0 - nf0) * rho_ss + nf0 * sl_max * rho_l_mid + nf0 * (1.0 - sl_max) * rho_g_mid;
        let sigma_v_mid = (hh - z) * g * rho_mid;
        if !with_gas {
            let tol = if incompressible_liq { 1e-11 } else { 1e-5 };
            assert_approx_eq!(rho_mid, (gg_s + e0) / (1.0 + e0), tol);
        }

        // println!("     pl_mid = {:>21} rho_l_mid = {:>21}", pl_mid, rho_l_mid);
        // println!("     pg_mid = {:>21} rho_g_mid = {:>21}", pg_mid, rho_g_mid);
        // println!("sigma_v_mid = {:>21}   rho_mid = {:>21}", sigma_v_mid, rho_mid);
        // println!();

        // parameters
        let footing = SampleParam::param_solid();
        let (upper, lower) = if with_gas {
            (
                SampleParam::param_porous_sol_liq_gas(nf0, 1e-2),
                SampleParam::param_porous_sol_liq_gas(0.2, 1e-2),
            )
        } else {
            (
                SampleParam::param_porous_sol_liq(nf0, 1e-2),
                SampleParam::param_porous_sol_liq(0.2, 1e-2),
            )
        };
        (footing, upper, lower, sigma_v_mid)
    }

    #[test]
    fn geostatics_new_works() -> Result<(), StrError> {
        // solid-liquid-gas
        let (footing, upper, lower, sigma_v_mid) = two_layers(3.0, 1.0, true, false);
        let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        config
            .elements(111, ElementConfig::Porous(lower, None))?
            .elements(222, ElementConfig::Porous(upper, None))?
            .elements(333, ElementConfig::Solid(footing, None))?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.layers.len(), 2);
        let top = &geo.layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.overburden, 0.0);
        let bottom = &geo.layers[1];
        assert_eq!(bottom.z_min, 0.0);
        assert_eq!(bottom.z_max, 1.0);
        assert_approx_eq!(bottom.overburden, sigma_v_mid, 1e-15);

        // solid-liquid
        let (footing, upper, lower, sigma_v_mid) = two_layers(3.0, 1.0, false, false);
        let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        config
            .elements(1, ElementConfig::Porous(lower, None))?
            .elements(2, ElementConfig::Porous(upper, None))?
            .elements(3, ElementConfig::Solid(footing, None))?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.layers.len(), 2);
        let top = &geo.layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.overburden, 0.0);
        let bottom = &geo.layers[1];
        assert_eq!(bottom.z_min, 0.0);
        assert_eq!(bottom.z_max, 1.0);
        assert_approx_eq!(bottom.overburden, sigma_v_mid, 1e-15);

        // solid-liquid(incompressible)
        let (footing, upper, lower, sigma_v_mid) = two_layers(3.0, 1.0, false, true);
        let mesh = Mesh::from_text_file("./data/meshes/column_two_layers_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        config
            .elements(1, ElementConfig::Porous(lower, None))?
            .elements(2, ElementConfig::Porous(upper, None))?
            .elements(3, ElementConfig::Solid(footing, None))?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.layers.len(), 2);
        let top = &geo.layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.overburden, 0.0);
        let bottom = &geo.layers[1];
        assert_eq!(bottom.z_min, 0.0);
        assert_eq!(bottom.z_max, 1.0);
        assert_approx_eq!(bottom.overburden, sigma_v_mid, 1e-15);
        Ok(())
    }
}
