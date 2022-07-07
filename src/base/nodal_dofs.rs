#![allow(unused)]

use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use gemlab::shapes::{GeoClass, GeoKind};
use std::collections::HashMap;
use std::fmt;

pub struct NodalDofs {
    pub all: HashMap<(Element, GeoKind), Vec<Vec<Dof>>>,
}

impl NodalDofs {
    pub fn new() -> Self {
        NodalDofs { all: HashMap::new() }
    }

    pub fn set(&mut self, ndim: usize, element: Element, kind: GeoKind) {
        // check consistency
        let rod_or_beam = element == Element::Rod || element == Element::Beam;
        let lin_geometry = kind.is_lin();
        if (rod_or_beam && !lin_geometry) || (!rod_or_beam && lin_geometry) {
            return; // inconsistent combination
        }

        // check for existent entry
        let key = &(element, kind);
        if self.all.contains_key(key) {
            return; // already configured
        }

        // update HashMap
        let nnode = kind.nnode();
        let dofs = match element {
            Element::Rod => {
                let dofs_per_node = if ndim == 2 {
                    vec![Dof::Ux, Dof::Uy]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz]
                };
                vec![dofs_per_node; nnode]
            }
            Element::Beam => {
                let dofs_per_node = if ndim == 2 {
                    vec![Dof::Ux, Dof::Uy, Dof::Rz]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
                };
                vec![dofs_per_node; nnode]
            }
            Element::Solid => {
                let dofs_per_node = if ndim == 2 {
                    vec![Dof::Ux, Dof::Uy]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz]
                };
                vec![dofs_per_node; nnode]
            }
            Element::PorousLiq => {
                let dofs_per_node = vec![Dof::Pl];
                vec![dofs_per_node; nnode]
            }
            Element::PorousLiqGas => {
                let dofs_per_node = vec![Dof::Pl, Dof::Pg];
                vec![dofs_per_node; nnode]
            }
            Element::PorousSldLiq => {
                let index = match POROUS_SLD_GEO_KIND_ALLOWED.iter().position(|k| *k == kind) {
                    Some(i) => i,
                    None => return, // not allowed
                };
                let high_order_dofs_per_node = if ndim == 2 {
                    vec![Dof::Ux, Dof::Uy]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz]
                };
                let mut dofs = vec![high_order_dofs_per_node; nnode];
                let low_order = kind.lower_order().unwrap();
                for m in 0..low_order.nnode() {
                    dofs[m].push(Dof::Pl);
                }
                dofs
            }
            Element::PorousSldLiqGas => {
                let index = match POROUS_SLD_GEO_KIND_ALLOWED.iter().position(|k| *k == kind) {
                    Some(i) => i,
                    None => return, // not allowed
                };
                let high_order_dofs_per_node = if ndim == 2 {
                    vec![Dof::Ux, Dof::Uy]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz]
                };
                let mut dofs = vec![high_order_dofs_per_node; nnode];
                let low_order = kind.lower_order().unwrap();
                for m in 0..low_order.nnode() {
                    dofs[m].push(Dof::Pl);
                    dofs[m].push(Dof::Pg);
                }
                dofs
            }
        };
        self.all.insert(*key, dofs);
    }

    pub fn validate(&self) -> Option<String> {
        for ((element, kind), dofs) in &self.all {
            //
        }
        None // all good
    }
}

impl fmt::Display for NodalDofs {
    /// Prints a formatted summary of Nodal Dofs
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Nodal Dofs\n").unwrap();
        write!(f, "==========\n").unwrap();
        let mut keys: Vec<_> = self.all.keys().collect();
        keys.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        for key in keys {
            let dofs = self.all.get(key).unwrap();
            write!(f, "{:?} → {:?}\n", key, dofs).unwrap();
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::NodalDofs;
    use crate::base::Element;
    use gemlab::shapes::GeoKind;

    #[test]
    pub fn nodal_dofs_work() {
        let mut nodal_dofs = NodalDofs::new();
        nodal_dofs.set(2, Element::Rod, GeoKind::Lin2);
        nodal_dofs.set(2, Element::Beam, GeoKind::Lin2);
        nodal_dofs.set(2, Element::Solid, GeoKind::Tri3);
        nodal_dofs.set(2, Element::PorousLiq, GeoKind::Tri3);
        nodal_dofs.set(2, Element::PorousLiqGas, GeoKind::Tri3);
        nodal_dofs.set(2, Element::PorousSldLiq, GeoKind::Tri6);
        nodal_dofs.set(2, Element::PorousSldLiqGas, GeoKind::Tri6);
        assert_eq!(format!("{}", nodal_dofs),
            "Nodal Dofs\n\
             ==========\n\
             (Rod, Lin2) → [[Ux, Uy], [Ux, Uy]]\n\
             (Beam, Lin2) → [[Ux, Uy, Rz], [Ux, Uy, Rz]]\n\
             (Solid, Tri3) → [[Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             (PorousLiq, Tri3) → [[Pl], [Pl], [Pl]]\n\
             (PorousLiqGas, Tri3) → [[Pl, Pg], [Pl, Pg], [Pl, Pg]]\n\
             (PorousSldLiq, Tri6) → [[Ux, Uy, Pl], [Ux, Uy, Pl], [Ux, Uy, Pl], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n\
             (PorousSldLiqGas, Tri6) → [[Ux, Uy, Pl, Pg], [Ux, Uy, Pl, Pg], [Ux, Uy, Pl, Pg], [Ux, Uy], [Ux, Uy], [Ux, Uy]]\n"
        );
    }
}
