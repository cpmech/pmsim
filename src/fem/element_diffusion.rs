use super::{Data, LocalEquations, State};
use crate::base::{Config, ParamDiffusion};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{copy_tensor2, Tensor2};

pub struct ElementDiffusion<'a> {
    pub data: &'a Data<'a>,
    pub config: &'a Config,
    pub cell: &'a Cell,
    pub param: &'a ParamDiffusion,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Matrix,
    pub conductivity: Tensor2,
}

impl<'a> ElementDiffusion<'a> {
    pub fn new(
        data: &'a Data,
        config: &'a Config,
        cell: &'a Cell,
        param: &'a ParamDiffusion,
    ) -> Result<Self, StrError> {
        // constants
        let ndim = data.mesh.ndim;
        let neq = data.element_dofs_map.get(cell)?.n_equation_local;

        // pad and ips
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(ndim, kind)?;
        set_pad_coords(&mut pad, &points, data.mesh);

        // conductivity
        let mut conductivity = Tensor2::new(true, ndim == 2);
        conductivity.sym_set(0, 0, param.kx);
        conductivity.sym_set(1, 1, param.ky);
        if ndim == 3 {
            conductivity.sym_set(2, 2, param.kz);
        }

        // done
        Ok({
            ElementDiffusion {
                data,
                config,
                cell,
                param,
                pad,
                ips: config.integ_point_data(cell)?,
                residual: Vector::new(neq),
                jacobian: Matrix::new(neq, neq),
                conductivity,
            }
        })
    }
}

impl<'a> LocalEquations for ElementDiffusion<'a> {
    fn residual(&mut self, state: &State) -> Result<(), StrError> {
        integ::vec_03_vg(&mut self.residual, &mut self.pad, 0, true, self.ips, |_w, _| {
            // todo
            Ok(())
        })?;
        if let Some(s) = self.param.source {
            integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, false, self.ips, |_| Ok(-s))?;
        }
        if self.config.control.transient {
            let theta = self.config.control.theta;
            let dt = state.delta_time;
            let alpha_1 = 1.0 / (theta * dt);
            integ::vec_01_ns(&mut self.residual, &mut self.pad, 0, false, self.ips, |_| {
                // TODO
                Ok(self.param.rho * (alpha_1 * 0.0))
            })?;
        }
        Ok(())
    }

    fn jacobian(&mut self, _state: &State) -> Result<(), StrError> {
        integ::mat_03_gtg(&mut self.jacobian, &mut self.pad, 0, 0, true, self.ips, |k, _| {
            copy_tensor2(k, &self.conductivity)
        })
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDiffusion;
    use crate::base::{Config, Element, ParamDiffusion};
    use crate::fem::{Data, LocalEquations, State};
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_chk::assert_vec_approx_eq;

    /*
    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let mut mesh_wrong = mesh.clone();
        mesh_wrong.cells[0].attribute_id = 100; // << never do this!

        let p1 = SampleParams::param_diffusion();
        let data = Data::new(&mesh_wrong, [(1, Element::Diffusion(p1))]).unwrap();
        let mut config = Config::new();
        assert_eq!(
            ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).err(),
            Some("cannot extract CellAttributeId to allocate ElementDiffusion")
        );

        config.n_integ_point.insert(1, 100); // wrong
        assert_eq!(
            ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).err(),
            Some("desired number of integration points is not available for Tri class")
        );
    }
    */

    // #[test]
    fn _element_diffusion_works() {
        let mesh = Samples::one_tri3();
        let rho = 1.0;
        let kx = 2.0;
        let ky = 3.0;
        let source = 4.0;
        let p1 = ParamDiffusion {
            rho,
            kx,
            ky,
            kz: 0.0,
            source: Some(source),
        };
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).unwrap();

        // check residual vector
        let mut state = State::new(&data, &config).unwrap();
        state.primary_unknowns[0] = 0.1;
        state.primary_unknowns[1] = 0.2;
        state.primary_unknowns[2] = 0.3;
        elem.residual(&state).unwrap();
        let ana = integ::AnalyticalTri3::new(&elem.pad);
        let (w0, w1) = (1.0, 2.0);
        let correct_r1 = ana.vec_01_ns(-source);
        let correct_r2 = ana.vec_03_vg(-w0, -w1);
        let correct_r = vec![
            correct_r1[0] + correct_r2[0],
            correct_r1[1] + correct_r2[1],
            correct_r1[2] + correct_r2[2],
        ];
        assert_vec_approx_eq!(elem.residual.as_data(), correct_r, 1e-15);

        // check Jacobian matrix
        elem.jacobian(&state).unwrap();
        let correct_kk = ana.mat_03_gtg(kx, ky);
        assert_vec_approx_eq!(elem.jacobian.as_data(), correct_kk.as_data(), 1e-12);
    }
}
