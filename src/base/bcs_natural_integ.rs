use super::{BcsNatural, FnBc, Nbc};
use crate::StrError;
use gemlab::integ;
use gemlab::mesh::{set_pad_coords, Mesh};
use gemlab::shapes::Scratchpad;
use russell_lab::Vector;

/// Returns the Scratchpads to perform the numerical integrations of the natural boundary conditions
pub fn get_pads(bcs: &BcsNatural, mesh: &Mesh) -> Result<Vec<(Scratchpad, Nbc, FnBc)>, StrError> {
    let mut results = Vec::new();
    for (face, nbc, f) in &bcs.faces {
        let mut pad = Scratchpad::new(mesh.ndim, face.kind)?;
        set_pad_coords(&mut pad, &face.points, &mesh);
        results.push((pad, *nbc, *f));
    }
    for (edge, nbc, f) in &bcs.edges {
        let mut pad = Scratchpad::new(mesh.ndim, edge.kind)?;
        set_pad_coords(&mut pad, &edge.points, &mesh);
        results.push((pad, *nbc, *f));
    }
    Ok(results)
}

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
