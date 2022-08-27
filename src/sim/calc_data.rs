#![allow(unused)]

use crate::base::{Config, DofNumbers, Element};
use crate::element::{self, Calculate};
use crate::StrError;
use gemlab::integ;
use gemlab::integ::{default_points, IntegPointData};
use gemlab::mesh::{set_pad_coords, Cell, CellId, Feature, Mesh};
use gemlab::shapes::{GeoKind, Scratchpad};
use rayon::prelude::*;
use russell_lab::{Matrix, Vector};

// pub type FnResidual = fn(residual: &mut Vector, time: f64, thickness: f64) -> Result<(), StrError>;
// pub type FnJacobian = fn(jacobian: &mut Matrix, time: f64, thickness: f64) -> Result<(), StrError>;

pub struct CalcData {
    // pub element: Element,
    pub pad: Scratchpad,
    pub ips: IntegPointData,
    pub residual: Vector,
    pub jacobian: Matrix,
    pub local_to_global: Vec<usize>,
    // pub fn_residual: FnResidual,
    // pub fn_jacobian: FnJacobian,
    // pub calculate: Box<dyn Calculate>,
}

impl CalcData {
    // Allocates new instance
    pub fn new(
        mesh: &Mesh,
        dn: &DofNumbers,
        config: &Config,
        cell_id: CellId,
        // element: Element,
        // calculate: Box<dyn Calculate>,
    ) -> Result<Self, StrError> {
        // pad and ips
        let cell = &mesh.cells[cell_id];
        let (kind, points) = (cell.kind, &cell.points);
        let mut pad = Scratchpad::new(mesh.ndim, kind).unwrap();
        set_pad_coords(&mut pad, &points, &mesh);
        let ips = default_points(pad.kind);

        // local_to_global
        let info = dn.element_dofs.get(&(cell.attribute_id, cell.kind)).unwrap();
        let mut local_to_global = vec![0; info.n_equation_local];
        for m in 0..points.len() {
            for (dof, local) in &info.dof_equation_pairs[m] {
                let global = *dn.point_dofs[cell.points[m]]
                    .get(dof)
                    .ok_or("cannot find DOF for CalcData")?;
                local_to_global[*local] = global;
            }
        }

        let element = dn.elements.get(&cell.attribute_id).unwrap();
        match element {
            Element::Rod(..) => (),
            Element::Beam(..) => (),
            Element::Solid(..) => (),
            Element::PorousLiq(..) => (),
            Element::PorousLiqGas(..) => (),
            Element::PorousSldLiq(..) => (),
            Element::PorousSldLiqGas(..) => (),
        }

        /*
        // function to compute the residual vector
        let fn_residual = match element {
            Element::Rod => element::Rod::fn_residual,
            Element::Beam => panic!("TODO"),
            Element::Solid => panic!("TODO"),
            Element::PorousLiq => panic!("TODO"),
            Element::PorousLiqGas => panic!("TODO"),
            Element::PorousSldLiq => panic!("TODO"),
            Element::PorousSldLiqGas => panic!("TODO"),
        };

        // function to compute the Jacobian matrix
        let fn_jacobian = match element {
            Element::Rod => element::Rod::fn_jacobian,
            Element::Beam => panic!("TODO"),
            Element::Solid => panic!("TODO"),
            Element::PorousLiq => panic!("TODO"),
            Element::PorousLiqGas => panic!("TODO"),
            Element::PorousSldLiq => panic!("TODO"),
            Element::PorousSldLiqGas => panic!("TODO"),
        };
        */

        // new instance
        Ok(CalcData {
            // element,
            pad,
            ips,
            residual: Vector::new(info.n_equation_local),
            jacobian: Matrix::new(info.n_equation_local, info.n_equation_local),
            local_to_global,
            // calculate,
            // fn_residual,
            // fn_jacobian,
        })
    }

    /// Calculates the residual vector at given time
    pub fn calc_residual(&mut self, time: f64, thickness: f64) -> Result<(), StrError> {
        // (self.fn_residual)(&mut self.residual, time, thickness)
        // self.calculate.residual()
        Err("not yet")
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::CalcData;
    use crate::base::{Config, DofNumbers, Element, SampleParams};
    use crate::element;
    use gemlab::mesh::Samples;
    use rayon::prelude::*;
    use std::collections::HashMap;

    #[test]
    fn calc_residual_works() {
        let mesh = Samples::three_tri3();
        let p1 = SampleParams::param_solid();
        let dn = DofNumbers::new(&mesh, HashMap::from([(1, Element::Solid(p1))])).unwrap();
        // let mut data: Vec<_> = mesh
        //     .cells
        //     .iter()
        //     .map(|cell| CalcData::new(&mesh, &dn, &mesh.cells[0], Element::Solid).unwrap())
        //     .collect();
        // data.par_iter_mut().for_each(|data| {
        //     data.calc_residual(0.0, 1.0).unwrap();
        // })
        let p1 = SampleParams::param_solid();
        let calculate = Box::new(element::Solid::new(p1));
        let config = Config::new();
        let calc_data = CalcData::new(&mesh, &dn, &config, 0).unwrap();
    }
}
