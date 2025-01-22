use super::Idealization;
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
/// * `ideal` -- The geometry idealization (axisymmetric, plane-strain, plane-stress, none)
/// * `l2g` -- The local to global map
/// * `ksi` -- The coordinate of the integration point (ξᵖ)
/// * `pad` -- Scratchpad to calculate interpolation functions
#[rustfmt::skip]
pub(crate) fn calculate_strain(
    eps: &mut Tensor2,
    uu: &Vector,
    ideal: &Idealization,
    l2g: &[usize],
    ksi: &[f64],
    pad: &mut Scratchpad,
) -> Result<(), StrError> {
    let nnode = pad.kind.nnode();
    pad.calc_gradient(ksi)?;
    let gg = &pad.gradient;
    eps.clear();
    if ideal.two_dim {
        for m in 0..nnode {
            eps.sym_add(0, 0, 1.0,  uu[l2g[0+2*m]] * gg.get(m,0));
            eps.sym_add(1, 1, 1.0,  uu[l2g[1+2*m]] * gg.get(m,1));
            eps.sym_add(0, 1, 1.0, (uu[l2g[0+2*m]] * gg.get(m,1) + uu[l2g[1+2*m]] * gg.get(m,0))/2.0);
        }
        if ideal.axisymmetric {
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
    use crate::base::{compute_local_to_global, Attributes, Config, ElementDofsMap, Equations, Etype, ParamSolid};
    use crate::base::{
        elastic_solution_horizontal_displacement_field, elastic_solution_shear_displacement_field,
        elastic_solution_vertical_displacement_field, generate_horizontal_displacement_field,
        generate_shear_displacement_field, generate_vertical_displacement_field,
    };
    use gemlab::mesh::Samples;
    use gemlab::shapes::Scratchpad;
    use russell_lab::vec_approx_eq;
    use russell_tensor::Tensor2;

    #[test]
    fn calc_delta_eps_works() {
        // irrelevant parameters here
        let young = 10_000.0; // kPa
        let poisson = 0.2; // [-]

        // strain magnitude (either ε_xx, ε_yy, or ε_xy)
        const STRAIN: f64 = 4.56;

        // loop over meshes
        for mesh in &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
            Samples::one_hex8(),
        ] {
            // incremental displacement field
            // (equal total displacements because initial displacements are zero)
            let duu_h = generate_horizontal_displacement_field(&mesh, STRAIN);
            let duu_v = generate_vertical_displacement_field(&mesh, STRAIN);
            let duu_s = generate_shear_displacement_field(&mesh, STRAIN);

            // solution
            let ndim = mesh.ndim;
            let (strain_h, _) = elastic_solution_horizontal_displacement_field(young, poisson, ndim, STRAIN);
            let (strain_v, _) = elastic_solution_vertical_displacement_field(young, poisson, ndim, STRAIN);
            let (strain_s, _) = elastic_solution_shear_displacement_field(young, poisson, ndim, STRAIN);

            // check the first cell/element only
            let cell = &mesh.cells[0];

            // local-to-global map
            let p1 = ParamSolid::sample_linear_elastic();
            let att = Attributes::from([(1, Etype::Solid(p1))]);
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
            let mut de = Tensor2::new(config.ideal.mandel());

            // check increment of strains for all integration points
            for p in 0..ips.npoint() {
                let iota = ips.coords(p);
                // horizontal strain
                calculate_strain(&mut de, &duu_h, &config.ideal, &l2g, iota, &mut pad).unwrap();
                vec_approx_eq(de.vector(), strain_h.vector(), 1e-13);
                // vertical strain
                calculate_strain(&mut de, &duu_v, &config.ideal, &l2g, iota, &mut pad).unwrap();
                vec_approx_eq(de.vector(), strain_v.vector(), 1e-14);
                // shear strain
                calculate_strain(&mut de, &duu_s, &config.ideal, &l2g, iota, &mut pad).unwrap();
                vec_approx_eq(de.vector(), strain_s.vector(), 1e-14);
            }
        }
    }
}
