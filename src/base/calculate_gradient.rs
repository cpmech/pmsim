use crate::StrError;
use gemlab::shapes::Scratchpad;
use russell_lab::Vector;

/// Calculates the gradient of a scalar quantity such as φ, pl, and pg at a Gauss point from the global U vector
///
/// Note: this function will work with any scalar quantity associated with the `l2g` map.
///
/// # Input
///
/// * `grad_phi` -- The gradient vector
/// * `uu` -- The global vector
/// * `l2g` -- The local to global map
/// * `ksi` -- The coordinate of the integration point (ξᵖ)
/// * `pad` -- Scratchpad to calculate interpolation functions
///
/// # Output
///
/// * `grad` -- Will contain the gradient vector
/// * This function also returns the scalar quantity interpolated to the Gauss point
pub(crate) fn calculate_gradient(
    grad_phi: &mut Vector,
    uu: &Vector,
    l2g: &[usize],
    ksi: &[f64],
    pad: &mut Scratchpad,
) -> Result<f64, StrError> {
    // shape function and its gradient
    (pad.fn_interp)(&mut pad.interp, ksi); // N
    pad.calc_gradient(ksi)?; // B
    let nn = &pad.interp;
    let bb = &pad.gradient;

    // constants
    let (space_ndim, nnode) = pad.xxt.dims();

    // interpolate ϕ at integration point
    let mut phi = 0.0;
    for m in 0..nnode {
        phi += nn[m] * uu[l2g[m]];
    }

    // interpolate ∇ϕ at integration point
    for i in 0..space_ndim {
        grad_phi[i] = 0.0;
        for m in 0..nnode {
            grad_phi[i] += bb.get(m, i) * uu[l2g[m]];
        }
    }

    Ok(phi)
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::calculate_gradient;
    use crate::base::{generate_scalar_field_ax_plus_by, ParamDiffusion, Schema};
    use gemlab::integ::Gauss;
    use gemlab::mesh::Samples;
    use russell_lab::{approx_eq, vec_approx_eq, Vector};

    #[test]
    fn calculate_gradient_works() {
        const A: f64 = 3.0;
        const B: f64 = 4.0;

        // loop over meshes
        for mesh in &[
            Samples::one_qua4(),
            Samples::three_tri3(),
            Samples::ring_eight_qua8_rad1_thick1(),
            Samples::one_hex8(),
        ] {
            let uu = generate_scalar_field_ax_plus_by(&mesh, A, B);

            // check the first cell/element only
            let cell = &mesh.cells[0];

            // local-to-global map
            let p1 = ParamDiffusion::sample();
            let mut schema = Schema::new();
            schema.add_diffusion(1, p1).build(&mesh).unwrap();
            let l2g = schema.local_to_global(cell.id).unwrap();

            // pad for numerical integration
            let mut pad = mesh.get_pad(cell.id);

            // integration points
            let gauss = Gauss::new(cell.kind);

            // gradient vector
            let mut grad_phi = Vector::new(mesh.ndim);

            // solution
            let mut x = Vector::new(mesh.ndim);
            let correct_grad = if mesh.ndim == 2 {
                Vector::from(&[A, B])
            } else {
                Vector::from(&[A, B, 0.0])
            };

            // check increment of strains for all integration points
            for p in 0..gauss.npoint() {
                let iota = gauss.coords(p);
                pad.calc_coords(&mut x, iota).unwrap();
                let phi = calculate_gradient(&mut grad_phi, &uu, &l2g, iota, &mut pad).unwrap();
                // println!("x = {:?}, phi = {:?}, grad_phi = {:?}", x.as_data(), phi, grad_phi.as_data());
                approx_eq(phi, A * x[0] + B * x[1], 1e-14);
                vec_approx_eq(&grad_phi, &correct_grad, 1e-13);
            }
        }
    }
}
