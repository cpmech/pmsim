#![allow(unused)]

use crate::base::{BcNatural, Conditions, DofNumbers, ParamSolid, ParamStressStrain};
use crate::StrError;
use gemlab::integ::{self, mat_10_gdg, IntegPointData};
use gemlab::mesh::{set_pad_coords, CellId, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::{copy_matrix, Matrix, Vector};
use russell_tensor::LinElasticity;

pub struct Solid {
    pub r_local: Vector,
    pub kk_local: Matrix,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub lin_elastic: bool,
}

impl Solid {
    pub fn new(
        mesh: &Mesh,
        cell_id: CellId,
        param: &ParamSolid,
        dn: &DofNumbers,
        cd: &Conditions,
    ) -> Result<Self, StrError> {
        let cell = &mesh.cells[cell_id];
        let neq = dn.local_to_global[cell_id].len();
        let r_local = Vector::new(neq);
        let mut kk_local = Matrix::new(neq, neq);
        let mut pad = Scratchpad::new(mesh.ndim, cell.kind)?;
        set_pad_coords(&mut pad, &cell.points, &mesh);
        let ips = integ::default_points(pad.kind);
        let lin_elastic = match param.stress_strain {
            ParamStressStrain::LinearElastic { young, poisson } => {
                let le = LinElasticity::new(young, poisson, param.two_dim, param.plane_stress);
                integ::mat_10_gdg(&mut kk_local, &mut pad, 0, 0, true, ips, |dd, _| {
                    // let in_array = model.get_modulus().mat.as_data();
                    // let out_array = dd.mat.as_mut_data();
                    // for i in 0..in_array.len() {
                    //     out_array[i] = th * in_array[i];
                    // }
                    // copy_matrix(&mut dd.mat, &le.get_modulus().mat)
                    Ok(())
                })?;
                true
            }
            _ => false,
        };
        Ok(Solid {
            r_local,
            kk_local,
            pad,
            ips,
            lin_elastic,
        })
    }

    pub fn calc_r_local(&mut self, bc: &BcNatural) {
        // todo
    }

    pub fn calc_kk_local(&mut self) {}
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::Solid;
    use crate::base::{Conditions, DofNumbers, Element, ParamSolid, ParamStressStrain, SampleMeshes};
    use russell_chk::assert_vec_approx_eq;
    use russell_lab::Matrix;
    use std::collections::HashMap;

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::bhatti_example_1dot6_bracket();
        let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid)])).unwrap();
        let param = ParamSolid {
            density: 1.0,
            stress_strain: ParamStressStrain::LinearElastic {
                young: 10_000.0,
                poisson: 0.2,
            },
            two_dim: true,
            plane_stress: true,
            thickness: 0.25,
            n_integ_point: None,
        };

        let cd = Conditions::new();

        let elem = Solid::new(&mesh, 0, &param, &dn, &cd).unwrap();
        assert_eq!(elem.r_local.dim(), 6);
        assert_eq!(elem.kk_local.dims(), (6, 6));

        // compare against results from Bhatti's book
        #[rustfmt::skip]
        let kk_bhatti = Matrix::from( &[
            [  9.765625000000001e+02,  0.000000000000000e+00, -9.765625000000001e+02,  2.604166666666667e+02,  0.000000000000000e+00, -2.604166666666667e+02],
            [  0.000000000000000e+00,  3.906250000000000e+02,  5.208333333333334e+02, -3.906250000000000e+02, -5.208333333333334e+02,  0.000000000000000e+00],
            [ -9.765625000000001e+02,  5.208333333333334e+02,  1.671006944444445e+03, -7.812500000000000e+02, -6.944444444444445e+02,  2.604166666666667e+02],
            [  2.604166666666667e+02, -3.906250000000000e+02, -7.812500000000000e+02,  2.126736111111111e+03,  5.208333333333334e+02, -1.736111111111111e+03],
            [  0.000000000000000e+00, -5.208333333333334e+02, -6.944444444444445e+02,  5.208333333333334e+02,  6.944444444444445e+02,  0.000000000000000e+00],
            [ -2.604166666666667e+02,  0.000000000000000e+00,  2.604166666666667e+02, -1.736111111111111e+03,  0.000000000000000e+00,  1.736111111111111e+03],
        ]);
        assert_vec_approx_eq!(elem.kk_local.as_data(), kk_bhatti.as_data(), 1e-12);
    }
}
