#![allow(unused)]

use crate::base::ParamPorous;
use gemlab::mesh::{CellId, Mesh};
use russell_lab::{Matrix, Vector};

pub struct PorousSldLiq {
    pub fe: Vector,
    pub ke: Matrix,
}

impl PorousSldLiq {
    pub fn new(mesh: &Mesh, cell_id: CellId, param: &ParamPorous) -> Self {
        PorousSldLiq {
            fe: Vector::new(0),
            ke: Matrix::new(0, 0),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::{ParamPorous, SampleMeshes, SampleParams};

    use super::PorousSldLiq;

    #[test]
    fn new_works() {
        let mesh = SampleMeshes::one_qua8();
        let param = SampleParams::param_porous_sol_liq(0.4, 1.0);
        let elem = PorousSldLiq::new(&mesh, 0, &param);
        assert_eq!(elem.fe.dim(), 0);
        assert_eq!(elem.ke.dims(), (0, 0));
    }
}
