use super::Config;
use crate::StrError;
use gemlab::shapes::Scratchpad;
use russell_lab::Vector;
use russell_tensor::Tensor2;

/// Calculates strain (ε) or strain increment (Δε) from the global (U) or (ΔU) vectors
///
/// # Input
///
/// * `eps` -- The (delta) strain tensor
/// * `uu` -- The global (delta) displacement vector
/// * `config` -- Configuration data
/// * `l2g` -- The local to global map
/// * `ksi` -- The coordinate of the integration point (ξᵖ)
/// * `pad` -- Scratchpad to calculate interpolation functions
#[rustfmt::skip]
pub(crate) fn calculate_strain(
    eps: &mut Tensor2,
    uu: &Vector,
    config: &Config,
    l2g: &[usize],
    ksi: &[f64],
    pad: &mut Scratchpad,
) -> Result<(), StrError> {
    let nnode = pad.kind.nnode();
    pad.calc_gradient(ksi)?;
    let gg = &pad.gradient;
    eps.clear();
    if config.two_dim {
        for m in 0..nnode {
            eps.sym_add(0, 0, 1.0,  uu[l2g[0+2*m]] * gg.get(m,0));
            eps.sym_add(1, 1, 1.0,  uu[l2g[1+2*m]] * gg.get(m,1));
            eps.sym_add(0, 1, 1.0, (uu[l2g[0+2*m]] * gg.get(m,1) + uu[l2g[1+2*m]] * gg.get(m,0))/2.0);
        }
        if config.axisymmetric {
            // calculate radius
            (pad.fn_interp)(&mut pad.interp, ksi);
            let nn = &pad.interp;
            let mut r = 0.0; // radius @ x(ξᵖ)
            for m in 0..nnode {
                r += nn[m] * pad.xxt.get(0, m);
            }
            // compute out-of-plane strain increment component
            for m in 0..nnode {
                eps.sym_add(2, 2, 1.0, uu[l2g[0 + 2 * m]] * nn[m] / r);
            }
        }
    } else {
        for m in 0..nnode {
            eps.sym_add(0, 0, 1.0,  uu[l2g[0+3*m]] * gg.get(m,0));
            eps.sym_add(1, 1, 1.0,  uu[l2g[1+3*m]] * gg.get(m,1));
            eps.sym_add(2, 2, 1.0,  uu[l2g[2+3*m]] * gg.get(m,2));
            eps.sym_add(0, 1, 1.0, (uu[l2g[0+3*m]] * gg.get(m,1) + uu[l2g[1+3*m]] * gg.get(m,0))/2.0);
            eps.sym_add(1, 2, 1.0, (uu[l2g[1+3*m]] * gg.get(m,2) + uu[l2g[2+3*m]] * gg.get(m,1))/2.0);
            eps.sym_add(0, 2, 1.0, (uu[l2g[0+3*m]] * gg.get(m,2) + uu[l2g[2+3*m]] * gg.get(m,0))/2.0);
        }
    }
    Ok(())
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calculate_strain;
    use crate::base::{compute_local_to_global, Attributes, Config, Element, ElementDofsMap, Equations, SampleParams};
    use crate::base::{
        generate_horizontal_displacement_field, generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use gemlab::mesh::Samples;
    use gemlab::shapes::Scratchpad;
    use russell_lab::vec_approx_eq;
    use russell_tensor::{Tensor2, SQRT_2};

    #[test]
    fn calc_delta_eps_works() {
        // loop over meshes
        for mesh in &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
            Samples::one_hex8(),
        ] {
            // incremental displacement field
            // (equal total displacements because initial displacements are zero)
            let strain = 4.56;
            let duu_h = generate_horizontal_displacement_field(&mesh, strain);
            let duu_v = generate_vertical_displacement_field(&mesh, strain);
            let duu_s = generate_shear_displacement_field(&mesh, strain);

            // correct increments of strain
            let ndim = mesh.ndim;
            let solution_h = if ndim == 2 {
                vec![strain, 0.0, 0.0, 0.0]
            } else {
                vec![strain, 0.0, 0.0, 0.0, 0.0, 0.0]
            };
            let solution_v = if ndim == 2 {
                vec![0.0, strain, 0.0, 0.0]
            } else {
                vec![0.0, strain, 0.0, 0.0, 0.0, 0.0]
            };
            let solution_s = if ndim == 2 {
                vec![0.0, 0.0, 0.0, strain * SQRT_2 / 2.0]
            } else {
                vec![0.0, 0.0, 0.0, strain * SQRT_2 / 2.0, 0.0, 0.0]
            };

            // check the first cell/element only
            let cell = &mesh.cells[0];

            // local-to-global map
            let p1 = SampleParams::param_solid();
            let att = Attributes::from([(1, Element::Solid(p1))]);
            let emap = ElementDofsMap::new(&mesh, &att).unwrap();
            let eqs = Equations::new(&mesh, &emap).unwrap();
            let l2g = compute_local_to_global(&emap, &eqs, cell).unwrap();

            // configuration
            let config = Config::new(&mesh);

            // pad for numerical integration
            let (kind, points) = (cell.kind, &cell.points);
            let mut pad = Scratchpad::new(ndim, kind).unwrap();
            mesh.set_pad(&mut pad, &points);

            // integration points
            let ips = config.integ_point_data(cell).unwrap();

            // strain increment
            let mut de = Tensor2::new(config.mandel);

            // check increment of strains for all integration points
            for p in 0..ips.len() {
                let ksi = &ips[p];
                // horizontal strain
                calculate_strain(&mut de, &duu_h, &config, &l2g, ksi, &mut pad).unwrap();
                vec_approx_eq(&de.vector(), &solution_h, 1e-13);
                // vertical strain
                calculate_strain(&mut de, &duu_v, &config, &l2g, ksi, &mut pad).unwrap();
                vec_approx_eq(&de.vector(), &solution_v, 1e-14);
                // shear strain
                calculate_strain(&mut de, &duu_s, &config, &l2g, ksi, &mut pad).unwrap();
                vec_approx_eq(&de.vector(), &solution_s, 1e-14);
            }
        }
    }
}
