use super::{FnBc, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::Mesh;
use gemlab::shapes::Scratchpad;
use russell_lab::Vector;

pub fn calculate(
    t: f64,
    thickness: f64,
    mesh: &Mesh,
    pad: &mut Scratchpad,
    nbc: Nbc,
    f: FnBc,
    ips: integ::IntegPointData,
) -> Result<(), StrError> {
    if mesh.ndim == 3 && nbc == Nbc::Qn {
        return Err("Qn natural boundary condition is not available for 3D edge");
    }
    let mut b = Vector::new(0);
    match nbc {
        Nbc::Qn => integ::vec_02_nv_bry(&mut b, pad, 0, true, ips, |v, _, un| {
            for i in 0..mesh.ndim {
                v[i] = thickness * f(t) * un[i];
            }
            Ok(())
        })?,
        Nbc::Qx => {
            // todo
        }
        Nbc::Qy => {
            // todo
        }
        Nbc::Qz => {
            // todo
        }
    }
    Ok(())
}
