use super::{layer::Layer, layer_info::LayerInfo};
use crate::{ElementConfig, SimConfig, StrError};
use gemlab::mesh::{At, CellAttributeId, CellId};
use russell_tensor::Tensor2;
use std::collections::{HashMap, HashSet};

/// Implements functions to compute geostatic stresses
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
/// * The gas phase is optional. In this case the gas pressure and
///   the intrinsic density of gas are assumed to be zero (atmospheric).
pub struct Geostatics {
    top_down_layers: Vec<Layer>,
}

impl Geostatics {
    /// Allocate a new instance
    ///
    /// Also discovers the layers from the Mesh, computes layers
    /// limits and allocates models.
    pub fn new(config: &SimConfig) -> Result<Self, StrError> {
        // mesh and space_ndim
        let mesh = config.mesh;
        let two_dim = mesh.space_ndim == 2;
        if mesh.space_ndim < 2 || mesh.space_ndim > 3 {
            return Err("Geostatics requires space_ndim = 2 or 3");
        }

        // param for fluids
        let param_fluids = match &config.param_fluids {
            Some(p) => p,
            None => return Err("Geostatics requires the definition of parameters for fluids"),
        };

        // find cells near x_min
        let x_min = mesh.coords_min[0];
        let mut cells_near_x_min: HashSet<CellId> = HashSet::new();
        if mesh.space_ndim == 2 {
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
        let index_z = mesh.space_ndim - 1;
        let mut infos: HashMap<CellAttributeId, LayerInfo> = HashMap::new();
        for cell_id in &cells_near_x_min {
            let cell = &mesh.cells[*cell_id];
            let element_config = config.get_element_config(cell.attribute_id)?;
            match element_config {
                ElementConfig::Porous(..) => (),
                _ => continue, // skip other types such as Rod, Beam, Solid, Seepage
            };
            let shape = mesh.alloc_shape_cell(*cell_id)?;
            match infos.get_mut(&cell.attribute_id) {
                Some(info) => {
                    info.z_min = f64::min(info.z_min, shape.coords_min[index_z]);
                    info.z_max = f64::max(info.z_max, shape.coords_max[index_z]);
                }
                None => {
                    infos.insert(
                        cell.attribute_id,
                        LayerInfo {
                            attribute_id: cell.attribute_id,
                            z_min: shape.coords_min[index_z],
                            z_max: shape.coords_max[index_z],
                        },
                    );
                }
            };
        }

        // sort layer infos
        let mut infos: Vec<_> = infos.keys().map(|k| *infos.get(k).unwrap()).collect();
        infos.sort_by(|a, b| b.z_min.partial_cmp(&a.z_min).unwrap());
        if infos.len() < 1 {
            return Err("Geostatics requires at least one layer == CellAttributeId");
        }

        // total height of the whole set of layers (column)
        let column_z_max = infos.first().unwrap().z_max;
        let column_z_min = infos.last().unwrap().z_min;
        let height = column_z_max - column_z_min;
        if column_z_min != 0.0 {
            return Err("the min elevation (y in 2D or z in 3D) of the mesh must be equal to 0.0");
        }
        if height <= 0.0 {
            return Err("the height of the mesh must be greater than zero");
        }

        // allocate layers
        let mut sigma_z_total_over = 0.0; // the first layer at the top has zero overburden
        let mut top_down_layers: Vec<Layer> = Vec::new();
        for info in infos {
            let element_config = config.get_element_config(info.attribute_id)?;
            let layer = match &element_config {
                ElementConfig::Porous(param_porous, _) => Layer::new(
                    info.attribute_id,
                    info.z_min,
                    info.z_max,
                    sigma_z_total_over,
                    param_fluids,
                    param_porous,
                    height,
                    config.gravity,
                    two_dim,
                )?,
                _ => panic!("INTERNAL ERROR: porous element_config is missing"), // not supposed to happen
            };
            sigma_z_total_over += layer.calc_sigma_z_total(layer.z_min)?;
            top_down_layers.push(layer);
        }

        // done
        Ok(Geostatics { top_down_layers })
    }

    /// Finds the cell attribute id of the layer containing a given elevation
    pub fn find_attribute_id(&self, z: f64) -> Result<CellAttributeId, StrError> {
        for layer in &self.top_down_layers {
            if z >= layer.z_min && z <= layer.z_max {
                return Ok(layer.attribute_id);
            }
        }
        Err("elevation is outside the porous region limits")
    }

    /// Calculates the geostatic total vertical stress at an elevation within the layer
    ///
    /// **Important:** Returns values using the continuum mechanics sign convention
    ///                where compression is negative.
    pub fn calc_sigma_z_total(&self, z: f64) -> Result<f64, StrError> {
        let layer = self.find_layer(z)?;
        layer.calc_sigma_z_total(z)
    }

    /// Returns the total or effective stress tensor at given elevation
    ///
    /// **Important:** Returns values using the continuum mechanics sign convention
    ///                where compression is negative.
    pub fn calc_stress(&self, z: f64, total_stress: bool) -> Result<Tensor2, StrError> {
        let layer = self.find_layer(z)?;
        layer.calc_stress(z, total_stress)
    }

    /// Finds the layer where a given elevation is located within
    fn find_layer(&self, z: f64) -> Result<&Layer, StrError> {
        for layer in &self.top_down_layers {
            if z >= layer.z_min && z <= layer.z_max {
                return Ok(layer);
            }
        }
        Err("elevation is outside the porous region limits")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Geostatics;
    use crate::{
        ElementConfig, ParamFluids, ParamPorous, ParamRealDensity, ParamSolid, SampleParam, SimConfig, StrError,
    };
    use gemlab::mesh::Mesh;
    use russell_chk::assert_approx_eq;

    // Returns the parameters of a two-layer porous column with a footing at the top and
    // the total vertical stress at the "middle" (transition between upper and lower layers)
    //
    // Returns: (fluids, footing, upper, lower, sigma_v_mid_approx)
    fn two_layers(
        height: f64,
        z_middle: f64,
        with_gas: bool,
        incompressible_liq: bool,
    ) -> (ParamFluids, ParamSolid, ParamPorous, ParamPorous, f64) {
        // constants
        let hh = height; // height of porous column
        let z = z_middle; // at the middle of upper layer
        let g = 10.0; // gravity
        let nf0 = 0.4; // initial porosity of top layer
        let ns = 1.0 - nf0; // partial fraction of solids of top layer
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
        let pl_mid = pl_ref + rho_l_ref * f64::exp_m1(cc_l * (hh - z) * g) / cc_l;
        let rho_l_mid = rho_l_ref + cc_l * (pl_mid - pl_ref);
        assert_approx_eq!(pl_mid, (hh - z) * g * rho_l_ref, 1e-4);
        assert_approx_eq!(rho_l_mid, rho_l_ref, 1e-5);

        // pressure and density: gas
        let pg_mid = pg_ref + rho_g_ref * f64::exp_m1(cc_g * (hh - z) * g) / cc_g;
        let rho_g_mid = rho_g_ref + cc_g * (pg_mid - pg_ref);
        assert_approx_eq!(pg_mid, (hh - z) * g * rho_g_ref, 1e-5);
        assert_approx_eq!(rho_g_mid, rho_g_ref, 1e-6);

        // partial density of the mixture and total vertical stress at the middle
        let sl = sl_max;
        let sg = 1.0 - sl;
        let nl = sl * nf0;
        let ng = sg * nf0;
        let rho_mid = ns * rho_ss + nl * rho_l_mid + ng * rho_g_mid;
        let sigma_v_mid_approx = -(hh - z) * g * rho_mid;
        if !with_gas {
            let tol = if incompressible_liq { 1e-11 } else { 1e-5 };
            assert_approx_eq!(rho_mid, (gg_s + e0) / (1.0 + e0), tol);
        }

        // println!("     pl_mid = {:>21} rho_l_mid = {:>21}", pl_mid, rho_l_mid);
        // println!("     pg_mid = {:>21} rho_g_mid = {:>21}", pg_mid, rho_g_mid);
        // println!("sigma_v_mid = {:>21}   rho_mid = {:>21}", sigma_v_mid_approx, rho_mid);
        // println!();

        // parameters
        let fluids = ParamFluids {
            density_liquid: ParamRealDensity {
                cc: cc_l,
                p_ref: pl_ref,
                rho_ref: rho_l_ref,
                tt_ref: 25.0,
            },
            density_gas: if with_gas {
                Some(ParamRealDensity {
                    cc: cc_g,
                    p_ref: pg_ref,
                    rho_ref: rho_g_ref,
                    tt_ref: 25.0,
                })
            } else {
                None
            },
        };
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
        (fluids, footing, upper, lower, sigma_v_mid_approx)
    }

    #[test]
    fn solid_liquid_incompressible_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text_file("./data/meshes/column_two_layers_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let (fluids, footing, upper, lower, sigma_v_mid_approx) = two_layers(3.0, 1.0, false, true);
        config
            .elements(1, ElementConfig::Porous(lower, None))?
            .elements(2, ElementConfig::Porous(upper, None))?
            .elements(3, ElementConfig::Solid(footing, None))?
            .set_param_fluids(fluids)?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.top_down_layers.len(), 2);
        let top = &geo.top_down_layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.calc_sigma_z_total(top.z_max)?, 0.0);
        assert_approx_eq!(top.calc_sigma_z_total(top.z_min)?, sigma_v_mid_approx, 1e-10);
        let bot = &geo.top_down_layers[1];
        assert_eq!(bot.z_min, 0.0);
        assert_eq!(bot.z_max, 1.0);
        assert_approx_eq!(bot.calc_sigma_z_total(bot.z_max)?, sigma_v_mid_approx, 1e-10);
        Ok(())
    }

    #[test]
    fn solid_liquid_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text_file("./data/meshes/column_distorted_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let (fluids, footing, upper, lower, sigma_v_mid_approx) = two_layers(3.0, 1.0, false, false);
        config
            .elements(1, ElementConfig::Porous(lower, None))?
            .elements(2, ElementConfig::Porous(upper, None))?
            .elements(3, ElementConfig::Solid(footing, None))?
            .set_param_fluids(fluids)?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.top_down_layers.len(), 2);
        let top = &geo.top_down_layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.calc_sigma_z_total(top.z_max)?, 0.0);
        assert_approx_eq!(top.calc_sigma_z_total(top.z_min)?, sigma_v_mid_approx, 1e-4); // only 1e-4 because of liquid compressibility
        let bot = &geo.top_down_layers[1];
        assert_eq!(bot.z_min, 0.0);
        assert_eq!(bot.z_max, 1.0);
        assert_approx_eq!(bot.calc_sigma_z_total(bot.z_max)?, sigma_v_mid_approx, 1e-4);
        Ok(())
    }

    #[test]
    fn solid_liquid_gas_works() -> Result<(), StrError> {
        let mesh = Mesh::from_text_file("./data/meshes/rectangle_tris_quads.msh")?;
        let mut config = SimConfig::new(&mesh);
        let (fluids, footing, upper, lower, sigma_v_mid_approx) = two_layers(3.0, 1.0, true, false);
        config
            .elements(111, ElementConfig::Porous(lower, None))?
            .elements(222, ElementConfig::Porous(upper, None))?
            .elements(333, ElementConfig::Solid(footing, None))?
            .set_param_fluids(fluids)?
            .set_gravity(10.0)?; // m/s²
        let geo = Geostatics::new(&config)?;
        assert_eq!(geo.top_down_layers.len(), 2);
        let top = &geo.top_down_layers[0];
        assert_eq!(top.z_min, 1.0);
        assert_eq!(top.z_max, 3.0);
        assert_eq!(top.calc_sigma_z_total(top.z_max)?, 0.0);
        assert_approx_eq!(top.calc_sigma_z_total(top.z_min)?, sigma_v_mid_approx, 1e-4); // only 1e-4 because of liquid compressibility
        let bot = &geo.top_down_layers[1];
        assert_eq!(bot.z_min, 0.0);
        assert_eq!(bot.z_max, 1.0);
        assert_approx_eq!(bot.calc_sigma_z_total(bot.z_max)?, sigma_v_mid_approx, 1e-4);
        Ok(())
    }
}
