use super::{Data, LocalEquations, State};
use crate::base::{compute_local_to_global, Config, ParamDiffusion};
use crate::model::{allocate_conductivity_model, Conductivity};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{t2_dot_vec, Tensor2};

/// Implements the local Diffusion Element equations
pub struct ElementDiffusion<'a> {
    /// Number of space dimensions
    pub ndim: usize,

    /// Global configuration
    pub config: &'a Config,

    /// The cell corresponding to this element
    pub cell: &'a Cell,

    /// Material parameters
    pub param: &'a ParamDiffusion,

    /// Local-to-global mapping
    pub local_to_global: Vec<usize>,

    /// Temporary variables for numerical integration
    pub pad: Scratchpad,

    /// Integration point coordinates and weights
    pub ips: integ::IntegPointData,

    /// Conductivity model
    pub model: Box<dyn Conductivity>,

    /// (temporary) Conductivity tensor at a single integration point
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
        let ndim = data.mesh.ndim;
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, data.mesh);
        Ok({
            ElementDiffusion {
                ndim,
                config,
                cell,
                param,
                local_to_global: compute_local_to_global(&data.information, &data.equations, cell)?,
                pad,
                ips: config.integ_point_data(cell)?,
                model: allocate_conductivity_model(&param.conductivity, ndim == 2),
                conductivity: Tensor2::new(true, ndim == 2),
                grad_tt: Vector::new(ndim),
            }
        })
    }
}

impl<'a> LocalEquations for ElementDiffusion<'a> {
    /// Returns the local-to-global mapping
    fn local_to_global(&self) -> &Vec<usize> {
        &self.local_to_global
    }

    /// Calculates the residual vector
    fn calc_residual(&mut self, residual: &mut Vector, state: &State) -> Result<(), StrError> {
        let ndim = self.ndim;
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, self.ips);
        args.alpha = self.config.thickness;
        args.axisymmetric = self.config.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        integ::vec_03_vb(residual, &mut args, |w, _, nn, gg| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..npoint {
                tt += nn[m] * state.uu[l2g[m]];
            }
            // interpolate ∇T at integration point
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..npoint {
                    self.grad_tt[i] += gg[m][i] * state.uu[l2g[m]];
                }
            }
            // compute conductivity tensor at integration point
            self.model.tensor(&mut self.conductivity, tt)?;
            // w must be negative as in the residual, however, w := -k.∇T
            // so the double negative is necessary to obtain -w = -(-k.∇T) = k.∇T
            t2_dot_vec(w, 1.0, &self.conductivity, &self.grad_tt)
        })
        .unwrap();

        // flag updates: very important from here on
        args.clear = false;

        if self.config.transient {
            // calculate alphas
            let (alpha_1, _) = self.config.control.alphas_transient(state.dt)?;
            let s = match self.param.source {
                Some(val) => val,
                None => 0.0,
            };

            // transient and source terms
            integ::vec_01_ns(residual, &mut args, |_, nn| {
                // interpolate T and T★ to integration point
                let (mut tt, mut tt_star) = (0.0, 0.0);
                for m in 0..npoint {
                    tt += nn[m] * state.uu[l2g[m]];
                    tt_star += nn[m] * state.uu_star[l2g[m]];
                }
                Ok(self.param.rho * (alpha_1 * tt - tt_star) - s)
            })
            .unwrap();
        } else {
            // source term only (steady case)
            if let Some(s) = self.param.source {
                integ::vec_01_ns(residual, &mut args, |_, _| Ok(-s)).unwrap();
            }
        }
        Ok(())
    }

    /// Calculates the Jacobian matrix
    fn calc_jacobian(&mut self, jacobian: &mut Matrix, state: &State) -> Result<(), StrError> {
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let mut args = integ::CommonArgs::new(&mut self.pad, self.ips);
        args.alpha = self.config.thickness;
        args.axisymmetric = self.config.axisymmetric;

        // conductivity term (always present, so we calculate it first with clear=true)
        integ::mat_03_btb(jacobian, &mut args, |kk, _, nn, _| {
            // interpolate T at integration point
            let mut tt = 0.0;
            for m in 0..npoint {
                tt += nn[m] * state.uu[l2g[m]];
            }
            // compute conductivity tensor at integration point
            self.model.tensor(kk, tt)
        })
        .unwrap();

        // very important from here on
        args.clear = false;

        // diffusion (mass) matrix
        if self.config.transient {
            let (alpha_1, _) = self.config.control.alphas_transient(state.dt)?;
            integ::mat_01_nsn(jacobian, &mut args, |_, _, _| Ok(self.param.rho * alpha_1)).unwrap();
        }
        Ok(())
    }

    /// Updates secondary values such as stresses and internal values
    ///
    /// Note that state.uu, state.vv, and state.aa have been updated already
    fn update_secondary_values(&mut self, _state: &State, _duu: &Vector) -> Result<(), StrError> {
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
    use gemlab::mesh::{Cell, Samples};
    use gemlab::shapes::GeoKind;
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
        config.n_integ_point.insert(1, 3);
        let wrong_cell = Cell {
            id: 0,
            attribute_id: 2, // wrong
            kind: GeoKind::Tri3,
            points: vec![0, 1, 2],
        };
        assert_eq!(
            ElementDiffusion::new(&data, &config, &wrong_cell, &p1).err(),
            Some("cannot find (CellAttributeId, GeoKind) in ElementInfoMap")
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
        let (kx, ky) = (0.1, 0.2);
        let w0 = -kx * dtt_dx;
        let w1 = 0.0;
        let correct_r = Vector::from(&ana.vec_03_vb(-w0, -w1));
        vec_approx_eq(residual.as_data(), correct_r.as_data(), 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let correct_kk = ana.mat_03_btb(kx, ky, false);
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
        let correct_src = ana.vec_01_ns(-source, false);
        let mut correct_r_new = Vector::new(neq);
        add_vectors(&mut correct_r_new, 1.0, &correct_r, 1.0, &correct_src).unwrap();
        vec_approx_eq(residual.as_data(), correct_r_new.as_data(), 1e-15);

        // error in transient -----------------------------------------------

        let mut config = Config::new();
        config.transient = true;
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).unwrap();
        state.dt = 0.0; // wrong
        assert_eq!(
            elem.calc_residual(&mut residual, &state).err(),
            Some("Δt is smaller than the allowed minimum")
        );
        assert_eq!(
            elem.calc_jacobian(&mut jacobian, &state).err(),
            Some("Δt is smaller than the allowed minimum")
        );
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
        let (kx, ky, kz) = (0.1, 0.2, 0.3);
        let w0 = -kx * dtt_dx;
        let w1 = 0.0;
        let w2 = -kz * dtt_dz;
        let correct_r = Vector::from(&ana.vec_03_vb(-w0, -w1, -w2));
        vec_approx_eq(residual.as_data(), correct_r.as_data(), 1e-15);

        // check Jacobian matrix
        let mut jacobian = Matrix::new(neq, neq);
        elem.calc_jacobian(&mut jacobian, &state).unwrap();
        let conductivity =
            Tensor2::from_matrix(&[[kx, 0.0, 0.0], [0.0, ky, 0.0], [0.0, 0.0, kz]], true, false).unwrap();
        let correct_kk = ana.mat_03_btb(&conductivity);
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
