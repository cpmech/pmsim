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
        integ::vec_03_vg(residual, pad, 0, true, self.ips, |w, _, gg| {
            // interpolate ∇T to integration point
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..npoint {
                    self.grad_tt[i] += gg[m][i] * state.uu[l2g[m]];
                }
            }
            // w must be negative as in the residual, however, w := -k.∇T
            // so the double negative is necessary to obtain -w = -(-k.∇T) = k.∇T
            t2_dot_vec(w, 1.0, &self.conductivity, &self.grad_tt)
        })?;
        if self.config.transient {
            let theta = self.config.control.theta;
            let alpha_1 = 1.0 / (theta * state.dt);
            let alpha_2 = (1.0 - theta) / theta;
            let s = match self.param.source {
                Some(val) => val,
                None => 0.0,
            };
            integ::vec_01_ns(residual, pad, 0, false, self.ips, |_, nn| {
                // interpolate T and T★ to integration point
                let (mut tt, mut tt_star) = (0.0, 0.0);
                for m in 0..npoint {
                    tt += nn[m] * state.uu[m];
                    tt_star += nn[m] * (alpha_1 * state.uu_old[m] + alpha_2 * state.vv_old[m]);
                }
                Ok(self.param.rho * (alpha_1 * tt - tt_star) - s)
            })?;
        } else {
            if let Some(s) = self.param.source {
                integ::vec_01_ns(residual, pad, 0, false, self.ips, |_, _| Ok(-s))?;
            }
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &State) -> Result<(), StrError> {
        integ::mat_03_gtg(jacobian, &mut self.pad, 0, 0, true, self.ips, |k, _, _| {
            copy_tensor2(k, &self.conductivity)
        })?;
        if self.config.transient {
            let theta = self.config.control.theta;
            let alpha_1 = 1.0 / (theta * state.dt);
            integ::mat_01_nsn(jacobian, &mut self.pad, 0, 0, false, self.ips, |_, _, _| {
                Ok(self.param.rho * alpha_1)
            })?;
        }
        Ok(())
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
    use russell_tensor::Tensor2;

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
        state.uu[0] = tt_field(mesh.points[0].coords[0]);
        state.uu[1] = tt_field(mesh.points[1].coords[0]);
        state.uu[2] = tt_field(mesh.points[2].coords[0]);

        // analytical solver
        let ana = integ::AnalyticalTri3::new(&elem.pad);

        // check residual vector
        let neq = 3;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
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
        let correct_src = Vector::from(&ana.vec_01_ns(-source));
        let mut correct_r_new = Vector::new(neq);
        add_vectors(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        vec_approx_eq(residual.as_data(), correct_r_new.as_data(), 1e-15);
    }

    #[test]
    fn element_diffusion_works_3d() {
        // mesh and parameters
        let mesh = Samples::one_tet4();
        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).unwrap();

        // set heat flow from the top to bottom and right to left
        let mut state = State::new(&data, &config).unwrap();
        let tt_field = |x, z| 100.0 + 7.0 * x + 3.0 * z;
        state.uu[0] = tt_field(mesh.points[0].coords[0], mesh.points[0].coords[2]);
        state.uu[1] = tt_field(mesh.points[1].coords[0], mesh.points[1].coords[2]);
        state.uu[2] = tt_field(mesh.points[2].coords[0], mesh.points[2].coords[2]);
        state.uu[3] = tt_field(mesh.points[3].coords[0], mesh.points[3].coords[2]);

        // analytical solver
        let ana = integ::AnalyticalTet4::new(&elem.pad);

        // check residual vector
        let neq = 4;
        let mut residual = Vector::new(neq);
        elem.calc_residual(&mut residual, &state).unwrap();
        let (dtt_dx, dtt_dz) = (7.0, 3.0);
        let w0 = -p1.kx * dtt_dx;
        let w1 = 0.0;
        let w2 = -p1.kz * dtt_dz;
        let correct_r = Vector::from(&ana.vec_03_vg(-w0, -w1, -w2));
        vec_approx_eq(residual.as_data(), correct_r.as_data(), 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let conductivity =
            Tensor2::from_matrix(&[[p1.kx, 0.0, 0.0], [0.0, p1.ky, 0.0], [0.0, 0.0, p1.kz]], true, false).unwrap();
        let correct_kk = ana.mat_03_gtg(&conductivity);
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
        let correct_src = Vector::from(&ana.vec_01_ns(-source));
        let mut correct_r_new = Vector::new(neq);
        add_vectors(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        vec_approx_eq(residual.as_data(), correct_r_new.as_data(), 1e-15);
    }
}
