use super::{Data, LocalEquations, State};
use crate::base::{Config, ParamDiffusion};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{copy_tensor2, t2_dot_vec, Tensor2};

/// Implements the local Diffusion Element equations
pub struct ElementDiffusion<'a> {
    /// Number of space dimensions
    pub ndim: usize,

    /// Local-to-global mapping
    pub local_to_global: &'a Vec<usize>,

    /// Global configuration
    pub config: &'a Config,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamDiffusion,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub ips: integ::IntegPointData,

    /// Conductivity tensor
    pub conductivity: Tensor2,

    /// (temporary) Gradient of temperature at a single integration point
    ///
    /// ∇T @ ip
    pub grad_tt: Vector,
}

impl<'a> ElementDiffusion<'a> {
    /// Allocates new instance
    pub fn new(
        data: &'a Data,
        config: &'a Config,
        cell: &'a Cell,
        param: &'a ParamDiffusion,
    ) -> Result<Self, StrError> {
        // scratchpad for numerical integration
        let ndim = data.mesh.ndim;
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, data.mesh);

        // conductivity
        let mut conductivity = Tensor2::new(true, ndim == 2);
        conductivity.sym_set(0, 0, param.kx);
        conductivity.sym_set(1, 1, param.ky);
        if ndim == 3 {
            conductivity.sym_set(2, 2, param.kz);
        }

        Ok({
            ElementDiffusion {
                ndim,
                local_to_global: &data.dof_numbers.local_to_global[cell.id],
                config,
                cell,
                param,
                pad,
                ips: config.integ_point_data(cell)?,
                conductivity,
                grad_tt: Vector::new(ndim),
            }
        })
    }
}

impl<'a> LocalEquations for ElementDiffusion<'a> {
    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &State) -> Result<(), StrError> {
        let ndim = self.ndim;
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let pad = &mut self.pad;
        let tt = &state.primary_unknowns;
        integ::vec_03_vg(residual, pad, 0, true, self.ips, |w, _, gg| {
            // interpolate ∇T to ip
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..npoint {
                    self.grad_tt[i] += gg[m][i] * tt[l2g[m]];
                }
            }
            // w must be negative as in the residual, however, w := -k.∇T
            // so the double negative is necessary to obtain -w = -(-k.∇T)
            t2_dot_vec(w, -(-1.0), &self.conductivity, &self.grad_tt)
        })?;
        if let Some(s) = self.param.source {
            integ::vec_01_ns(residual, pad, 0, false, self.ips, |_, _| Ok(-s))?;
        }
        if self.config.control.transient {
            let theta = self.config.control.theta;
            let dt = state.delta_time;
            let alpha_1 = 1.0 / (theta * dt);
            integ::vec_01_ns(residual, pad, 0, false, self.ips, |_, _| {
                // TODO
                Ok(self.param.rho * (alpha_1 * 0.0))
            })?;
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, _state: &State) -> Result<(), StrError> {
        integ::mat_03_gtg(jacobian, &mut self.pad, 0, 0, true, self.ips, |k, _, _| {
            copy_tensor2(k, &self.conductivity)
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDiffusion;
    use crate::base::{Config, Element, SampleParams};
    use crate::fem::{Data, LocalEquations, State};
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_chk::vec_approx_eq;
    use russell_lab::{add_vectors, Matrix, Vector};

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let mut config = Config::new();
        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }

    #[test]
    fn element_diffusion_works_2d() {
        // mesh and parameters
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).unwrap();

        // set heat flow from the right to the left
        let mut state = State::new(&data, &config).unwrap();
        let tt_field = |x| 100.0 + 5.0 * x;
        state.primary_unknowns[0] = tt_field(mesh.points[0].coords[0]);
        state.primary_unknowns[1] = tt_field(mesh.points[1].coords[0]);
        state.primary_unknowns[2] = tt_field(mesh.points[2].coords[0]);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
        let neq = 3;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        println!("{}", residual);
        let dtt_dx = 5.0;
        let w0 = -p1.kx * dtt_dx;
        let w1 = 0.0;
        let correct_r = Vector::from(&ana.vec_03_vg(-w0, -w1));
        vec_approx_eq(residual.as_data(), correct_r.as_data(), 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct_kk = ana.mat_03_gtg(p1.kx, p1.ky);
        vec_approx_eq(jacobian.as_data(), correct_kk.as_data(), 1e-15);

        // with source term -------------------------------------------------

        // parameters
        let source = 4.0;
        let mut p1_new = p1.clone();
        p1_new.source = Some(source);
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1_new))]).unwrap();
        let config = Config::new();
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1_new).unwrap();

        // check residual vector
        elem.calc_residual(&mut residual, &state).unwrap();
        println!("{}", residual);
        let correct_src = Vector::from(&ana.vec_01_ns(-source));
        let mut correct_r_new = Vector::new(neq);
        add_vectors(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        println!("{}", correct_r_new);
        vec_approx_eq(residual.as_data(), correct_r_new.as_data(), 1e-15);
    }
}
