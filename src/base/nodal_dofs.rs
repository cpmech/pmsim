use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::shapes::GeoKind;
use std::collections::HashMap;
use std::fmt;

/// Holds the nodal degrees-of-freedom for a pair of (Element,GeoKind)
///
/// Example of DOFs for ([super::Element::PorousSldLiq],[gemlab::shapes::Tri6]):
///
/// ```text
///              {Ux,Uy,Pl}
///                  2
///                 / \
///                /   \
///       {Ux,Uy} 5     4 {Ux,Uy}
///              /       \
///             /         \
/// {Ux,Uy,Pl} 0-----3-----1 {Ux,Uy,Pl}
///               {Ux,Uy}
/// ```
pub struct NodalDofs {
    /// Space number of dimension (for verification purposes)
    ndim: usize,

    /// Holds all combinations
    pub all: HashMap<(Element, GeoKind), Vec<Vec<Dof>>>,
}

impl NodalDofs {
    /// Allocates a new instance
    pub fn new(ndim: usize) -> Self {
        NodalDofs {
            ndim,
            all: HashMap::new(),
        }
    }

    /// Sets the nodal DOFs for a pair of (Element,GeoKind)
    pub fn set(&mut self, element: Element, kind: GeoKind) -> Result<(), StrError> {
        // check consistency
        let rod_or_beam = element == Element::Rod || element == Element::Beam;
        let lin_geometry = kind.is_lin();
        if rod_or_beam && !lin_geometry {
            return Err("cannot set Rod or Beam with a non-Lin GeoClass"); // inconsistent combination
        }
        if !rod_or_beam && lin_geometry {
            return Err("GeoClass::Lin is reserved for Rod or Beam"); // inconsistent combination
        }

        // check for existent entry
        let key = &(element, kind);
        if self.all.contains_key(key) {
            return Ok(()); // already configured
        }

        // update HashMap
        let nnode = kind.nnode();
        let dofs = match element {
            Element::Rod => {
                let dofs_per_node = if self.ndim == 2 {
                    vec![Dof::Ux, Dof::Uy]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz]
                };
                vec![dofs_per_node; nnode]
            }
            Element::Beam => {
                let dofs_per_node = if self.ndim == 2 {
                    vec![Dof::Ux, Dof::Uy, Dof::Rz]
                } else {
                    vec![Dof::Ux, Dof::Uy, Dof::Uz, Dof::Rx, Dof::Ry, Dof::Rz]
                };
                vec![dofs_per_node; nnode]
            }
            Element::Solid => {
                let dofs_per_node = if self.ndim == 2 {
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
                if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                    return Err("cannot set PorousSldLiq with given GeoKind");
                };
                let high_order_dofs_per_node = if self.ndim == 2 {
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
                if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                    return Err("cannot set PorousSldLiqGas with given GeoKind");
                };
                let high_order_dofs_per_node = if self.ndim == 2 {
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
        Ok(())
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
    use crate::base::{Dof, Element};
    use crate::StrError;
    use gemlab::shapes::GeoKind;

    #[test]
    pub fn set_captures_errors() {
        let mut nodal_dofs = NodalDofs::new(2);
        assert_eq!(
            nodal_dofs.set(Element::Rod, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            nodal_dofs.set(Element::Beam, GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        assert_eq!(
            nodal_dofs.set(Element::Solid, GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        assert_eq!(
            nodal_dofs.set(Element::PorousSldLiq, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        assert_eq!(
            nodal_dofs.set(Element::PorousSldLiqGas, GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    pub fn set_and_get_work_2d() -> Result<(), StrError> {
        let mut nodal_dofs = NodalDofs::new(2);
        nodal_dofs.set(Element::Rod, GeoKind::Lin2)?;
        nodal_dofs.set(Element::Beam, GeoKind::Lin2)?;
        nodal_dofs.set(Element::Solid, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousLiq, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousLiq, GeoKind::Tri3)?; // already configured
        nodal_dofs.set(Element::PorousLiqGas, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousSldLiq, GeoKind::Tri6)?;
        nodal_dofs.set(Element::PorousSldLiqGas, GeoKind::Tri6)?;
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
        assert_eq!(
            nodal_dofs.all.get(&(Element::Rod, GeoKind::Lin2)),
            Some(&vec![vec![Dof::Ux, Dof::Uy], vec![Dof::Ux, Dof::Uy]])
        );
        assert_eq!(nodal_dofs.all.get(&(Element::PorousLiq, GeoKind::Qua4)), None);
        Ok(())
    }

    #[test]
    pub fn set_and_get_work_3d() -> Result<(), StrError> {
        let mut nodal_dofs = NodalDofs::new(3);
        nodal_dofs.set(Element::Rod, GeoKind::Lin2)?;
        nodal_dofs.set(Element::Beam, GeoKind::Lin2)?;
        nodal_dofs.set(Element::Solid, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousLiq, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousLiq, GeoKind::Tri3)?; // already configured
        nodal_dofs.set(Element::PorousLiqGas, GeoKind::Tri3)?;
        nodal_dofs.set(Element::PorousSldLiq, GeoKind::Tri6)?;
        nodal_dofs.set(Element::PorousSldLiqGas, GeoKind::Tri6)?;
        assert_eq!(format!("{}", nodal_dofs),
            "Nodal Dofs\n\
             ==========\n\
             (Rod, Lin2) → [[Ux, Uy, Uz], [Ux, Uy, Uz]]\n\
             (Beam, Lin2) → [[Ux, Uy, Uz, Rx, Ry, Rz], [Ux, Uy, Uz, Rx, Ry, Rz]]\n\
             (Solid, Tri3) → [[Ux, Uy, Uz], [Ux, Uy, Uz], [Ux, Uy, Uz]]\n\
             (PorousLiq, Tri3) → [[Pl], [Pl], [Pl]]\n\
             (PorousLiqGas, Tri3) → [[Pl, Pg], [Pl, Pg], [Pl, Pg]]\n\
             (PorousSldLiq, Tri6) → [[Ux, Uy, Uz, Pl], [Ux, Uy, Uz, Pl], [Ux, Uy, Uz, Pl], [Ux, Uy, Uz], [Ux, Uy, Uz], [Ux, Uy, Uz]]\n\
             (PorousSldLiqGas, Tri6) → [[Ux, Uy, Uz, Pl, Pg], [Ux, Uy, Uz, Pl, Pg], [Ux, Uy, Uz, Pl, Pg], [Ux, Uy, Uz], [Ux, Uy, Uz], [Ux, Uy, Uz]]\n"
        );
        assert_eq!(
            nodal_dofs.all.get(&(Element::Rod, GeoKind::Lin2)),
            Some(&vec![vec![Dof::Ux, Dof::Uy, Dof::Uz], vec![Dof::Ux, Dof::Uy, Dof::Uz]])
        );
        assert_eq!(nodal_dofs.all.get(&(Element::PorousLiq, GeoKind::Qua4)), None);
        Ok(())
    }
}
