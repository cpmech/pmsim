#![allow(unused)]

use super::{Data, ElementTrait, State};
use crate::base::{Config, ParamDiffusion};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Cell};
use gemlab::shapes::Scratchpad;
use russell_lab::{Matrix, Vector};
use russell_tensor::{copy_tensor2, t2_dot_vec, Tensor2};

pub struct ElementDiffusion<'a> {
    pub ndim: usize,
    pub local_to_global: &'a Vec<usize>,
    pub config: &'a Config,
    pub cell: &'a Cell,
    pub param: &'a ParamDiffusion,
    pub pad: Scratchpad,
    pub ips: integ::IntegPointData,
    pub residual: Vector,
    pub jacobian: Matrix,
    pub conductivity: Tensor2,

    pub grad_tt: Vector, // ∇T @ ip
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
        let neq = data.element_dofs.get(cell).unwrap().n_equation_local;

        // pad and ips
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

        // done
        Ok({
            ElementDiffusion {
                ndim,
                local_to_global: &data.dof_numbers.local_to_global[cell.id],
                config,
                cell,
                param,
                pad,
                ips: config.integ_point_data(cell)?,
                residual: Vector::new(neq),
                jacobian: Matrix::new(neq, neq),
                conductivity,
                grad_tt: Vector::new(ndim),
            }
        })
    }
}

impl<'a> ElementTrait for ElementDiffusion<'a> {
    fn residual(&mut self, state: &State) -> Result<(), StrError> {
        let ndim = self.ndim;
        let npoint = self.cell.points.len();
        let l2g = &self.local_to_global;
        let res = &mut self.residual;
        let pad = &mut self.pad;
        let tt = &state.primary_unknowns;
        integ::vec_03_vg(res, pad, 0, true, self.ips, |w, _, _, gg| {
            for i in 0..ndim {
                self.grad_tt[i] = 0.0;
                for m in 0..npoint {
                    self.grad_tt[i] += gg[m][i] * tt[l2g[m]];
                }
            }
            t2_dot_vec(w, -1.0, &self.conductivity, &self.grad_tt)
        })?;
        if let Some(s) = self.param.source {
            integ::vec_01_ns(res, pad, 0, false, self.ips, |_| Ok(-s))?;
        }
        if self.config.control.transient {
            let theta = self.config.control.theta;
            let dt = state.delta_time;
            let alpha_1 = 1.0 / (theta * dt);
            integ::vec_01_ns(res, pad, 0, false, self.ips, |_| {
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
    use crate::base::{Config, Element, ParamDiffusion, SampleMeshes, SampleParams};
    use crate::fem::{Data, ElementTrait, State};
    use gemlab::integ;
    use gemlab::mesh::Samples;
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::{copy_vector, Matrix, Vector};

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
    fn element_diffusion_works_bhatti_6dot22() {
        let mesh = SampleMeshes::bhatti_example_6dot22_heat();
        let kx = 45.0;
        let ky = 45.0;
        let source = 5e6;
        let p1 = ParamDiffusion {
            rho: 0.0,
            kx,
            ky,
            kz: 0.0,
            source: Some(source),
        };
        let data = Data::new(&mesh, [(1, Element::Diffusion(p1))]).unwrap();
        let config = Config::new();
        let mut elem = ElementDiffusion::new(&data, &config, &mesh.cells[0], &p1).unwrap();

        // check residual of first element
        // with state = 0, the residual is equal to -b (negative of integral of source term)
        let mut state = State::new(&data, &config).unwrap();
        let neg_b = Vector::from(&[250.0, 312.5, 312.5, 250.0, -1125., -1000.0, -1125., -1250.0]);
        elem.residual(&state).unwrap();
        assert_vec_approx_eq!(elem.residual.as_data(), neg_b.as_data(), 1e-12);

        // check Jacobian of first element (independent of state)
        elem.jacobian(&state).unwrap();
        #[rustfmt::skip]
        let bhatti_kk = Matrix::from(&[
            [38.515873015873005  , 23.194444444444443  , 21.46825396825396   , 22.02777777777777   , -20.317460317460306 , -30.91269841269842  , -14.682539682539685 , -39.293650793650784],
            [23.19444444444444   , 84.08730158730152   , 33.6111111111111    , 30.456349206349195  , -26.984126984126974 , -82.69841269841265  , -28.015873015873005 , -33.650793650793645],
            [21.46825396825396   , 33.6111111111111    , 58.313492063492056  , 23.19444444444444   , -19.36507936507935  , -68.17460317460318  , -20.63492063492062  , -28.412698412698408],
            [22.02777777777777   , 30.456349206349195  , 23.19444444444444   , 57.14682539682538   , -16.507936507936492 , -36.150793650793645 , -33.49206349206349  , -46.67460317460315 ],
            [-20.317460317460302 , -26.984126984126966 , -19.365079365079353 , -16.507936507936492 , 95.07936507936506   , 16.349206349206334  , -5.0793650793650915 , -23.174603174603188],
            [-30.91269841269841  , -82.69841269841265  , -68.17460317460318  , -36.15079365079365  , 16.34920634920633   , 155.87301587301585  , 3.650793650793659   , 42.06349206349206  ],
            [-14.682539682539685 , -28.015873015873005 , -20.63492063492062  , -33.49206349206349  , -5.079365079365089  , 3.650793650793652   , 95.07936507936506   , 3.174603174603193  ],
            [-39.293650793650784 , -33.65079365079365  , -28.412698412698408 , -46.67460317460315  , -23.174603174603188 , 42.06349206349206   , 3.1746031746031935  , 125.96825396825392 ],
        ]);
        assert_vec_approx_eq!(elem.jacobian.as_data(), bhatti_kk.as_data(), 1e-13);

        return;

        let tt_bhatti = Vector::from(&[
            156.4451416316945,
            150.75534475609425,
            149.20088556384715,
            144.22821336874946,
            133.84330164705196,
            123.99835257095894,
            121.7541135843015,
            119.15086573652133,
            110.0,
            110.0,
            110.0,
            144.67722333158,
            129.13353090328047,
        ]);
        copy_vector(&mut state.primary_unknowns, &tt_bhatti);
        elem.residual(&state).unwrap();
        println!("{}", elem.residual);
        // println!("{}", elem.jacobian);
    }

    #[test]
    fn element_diffusion_works() {
        let mesh = Samples::one_tri3();
        let rho = 1.0;
        let kx = 2.0;
        let ky = 3.0;
        let p1 = ParamDiffusion {
            rho,
            kx,
            ky,
            kz: 0.0,
            source: None,
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
        let correct_r = ana.vec_03_vg(-w0, -w1);
        println!("{}", elem.residual);
        println!("{:?}", correct_r);
        // assert_vec_approx_eq!(elem.residual.as_data(), correct_r, 1e-15);

        // check Jacobian matrix
        elem.jacobian(&state).unwrap();
        let correct_kk = ana.mat_03_gtg(kx, ky);
        // assert_vec_approx_eq!(elem.jacobian.as_data(), correct_kk.as_data(), 1e-12);

        // let source = 4.0;
        // let p1 = ParamDiffusion {
        //     rho,
        //     kx,
        //     ky,
        //     kz: 0.0,
        //     source: Some(source),
        // };
        // let correct_r1 = ana.vec_01_ns(-source);
        // let correct_r = vec![
        //     correct_r1[0] + correct_r2[0],
        //     correct_r1[1] + correct_r2[1],
        //     correct_r1[2] + correct_r2[2],
        // ];
    }
}
