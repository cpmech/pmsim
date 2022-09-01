use super::{Dof, Element, POROUS_SLD_GEO_KIND_ALLOWED};
use crate::StrError;
use gemlab::mesh::{CellAttributeId, Mesh};
use gemlab::shapes::GeoKind;
use std::collections::HashMap;
use std::fmt::{self, Write};

/// Holds the DOFs and local equation numbers of an Element/GeoKind pair
///
/// ```text
/// leq: local equation number       leq   point   geq
/// geq: global equation number       ↓        ↓    ↓
///                                   0 → Ux @ 0 →  0
///            {Ux → 6}               1 → Uy @ 0 →  1
///            {Uy → 7}               2 → Ux @ 1 →  3
///            {Pl → 8}               3 → Uy @ 1 →  4
///                2                  4 → Ux @ 2 →  6
///               / \                 5 → Uy @ 2 →  7
///   {Ux → 13}  /   \  {Ux → 11}     6 → Ux @ 3 →  9
///   {Uy → 14} 5     4 {Uy → 12}     7 → Uy @ 3 → 10
///            /       \              8 → Ux @ 4 → 11
/// {Ux → 0}  /         \  {Ux → 3}   9 → Uy @ 4 → 12
/// {Uy → 1} 0-----3-----1 {Uy → 4}  10 → Ux @ 5 → 13
/// {Pl → 2}   {Ux → 9}    {Pl → 5}  11 → Uy @ 5 → 14
///            {Uy → 10}             12 → Pl @ 0 →  2  <<< eq_first_pl
///                                  13 → Pl @ 1 →  5
///                                  14 → Pl @ 2 →  8
/// ```
pub struct ElementDofs {
    /// Holds all cell DOF keys and local equation numbers
    ///
    /// **Notes:** The outer array has length = nnode.
    /// The inner arrays have variable lengths = ndof at the node.
    /// The inner arrays contain pairs of Dof and local equation numbers.
    pub dof_equation_pairs: Vec<Vec<(Dof, usize)>>,

    /// Dimension of the local system of equations
    ///
    /// **Note:** This is equal to the total number of DOFs in the cell
    pub n_equation_local: usize,

    /// Local equation number of the first Dof::Pl
    pub eq_first_pl: Option<usize>,

    /// Local equation number of the first Dof::Pg
    pub eq_first_pg: Option<usize>,

    /// Local equation number of the first Dof::T
    pub eq_first_tt: Option<usize>,
}

impl ElementDofs {
    /// Allocates a new instance
    pub fn new(ndim: usize, element: Element, kind: GeoKind) -> Result<Self, StrError> {
        // check
        let rod_or_beam = match element {
            Element::Rod(..) => true,
            Element::Beam(..) => true,
            _ => false,
        };
        let lin_geometry = kind.is_lin();
        if rod_or_beam && !lin_geometry {
            return Err("cannot set Rod or Beam with a non-Lin GeoClass"); // inconsistent combination
        }
        if !rod_or_beam && lin_geometry {
            return Err("GeoClass::Lin is reserved for Rod or Beam"); // inconsistent combination
        }

        // auxiliary data
        let nnode = kind.nnode();
        let mut dofs = vec![Vec::new(); nnode];
        let mut count = 0;
        let mut eq_first_pl = None;
        let mut eq_first_pg = None;
        let eq_first_tt = None;

        // handle each combination
        #[rustfmt::skip]
        match element {
            Element::Diffusion(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::T, count)); count += 1;
                }
            }
            Element::Rod(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
            }
            Element::Beam(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 2 {
                        dofs[m].push((Dof::Rz, count)); count += 1;
                    } else {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                        dofs[m].push((Dof::Rx, count)); count += 1;
                        dofs[m].push((Dof::Ry, count)); count += 1;
                        dofs[m].push((Dof::Rz, count)); count += 1;
                    }
                }
            }
            Element::Solid(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
            }
            Element::PorousLiq(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
            }
            Element::PorousLiqGas(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                    dofs[m].push((Dof::Pg, count)); count += 1;
                }
            }
            Element::PorousSldLiq(..) => {
                if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                    return Err("cannot set PorousSldLiq with given GeoKind");
                };
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
                let ncorner = kind.lower_order().unwrap().nnode();
                eq_first_pl = Some(count);
                for m in 0..ncorner {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
            }
            Element::PorousSldLiqGas(..) => {
                if !POROUS_SLD_GEO_KIND_ALLOWED.contains(&kind) {
                    return Err("cannot set PorousSldLiqGas with given GeoKind");
                };
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
                let ncorner = kind.lower_order().unwrap().nnode();
                eq_first_pl = Some(count);
                for m in 0..ncorner {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
                eq_first_pg = Some(count);
                for m in 0..ncorner {
                    dofs[m].push((Dof::Pg, count)); count += 1;
                }
            }
        };
        Ok(ElementDofs {
            dof_equation_pairs: dofs,
            n_equation_local: count,
            eq_first_pl,
            eq_first_pg,
            eq_first_tt,
        })
    }

    /// Allocates a new collection of ElementDofs
    pub fn new_collection(
        mesh: &Mesh,
        elements: &HashMap<CellAttributeId, Element>,
    ) -> Result<HashMap<(CellAttributeId, GeoKind), ElementDofs>, StrError> {
        let mut collection = HashMap::new();
        for cell in &mesh.cells {
            let element = match elements.get(&cell.attribute_id) {
                Some(e) => e,
                None => return Err("cannot find CellAttributeId in elements map to create new ElementDofs collection"),
            };
            collection.insert(
                (cell.attribute_id, cell.kind),
                ElementDofs::new(mesh.ndim, *element, cell.kind)?,
            );
        }
        Ok(collection)
    }
}

impl fmt::Display for ElementDofs {
    /// Displays an ElementDofs
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for m in 0..self.dof_equation_pairs.len() {
            write!(f, "{}: {:?}\n", m, self.dof_equation_pairs[m]).unwrap();
        }
        write!(
            f,
            "(Pl @ {:?}, Pg @ {:?}, T @ {:?})\n",
            self.eq_first_pl, self.eq_first_pg, self.eq_first_tt
        )
    }
}

/// Returns a string to display a collection of ElementDofs
pub fn string_element_dofs_collection(
    elements: &HashMap<CellAttributeId, Element>,
    collection: &HashMap<(CellAttributeId, GeoKind), ElementDofs>,
) -> String {
    let mut b = String::new();
    write!(&mut b, "Elements: DOFs and local equation numbers\n").unwrap();
    write!(&mut b, "=========================================\n").unwrap();
    let mut keys: Vec<_> = collection.keys().collect();
    keys.sort_by(|a, b| a.0.cmp(&b.0));
    for key in keys {
        let element_dofs = collection.get(key).unwrap();
        let (attr, kind) = key;
        let element = elements.get(attr).unwrap();
        write!(&mut b, "{} → {} → {:?}\n", attr, element.name(), kind).unwrap();
        write!(&mut b, "{}", element_dofs).unwrap();
        write!(&mut b, "-----------------------------------------\n").unwrap();
    }
    b
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use std::collections::HashMap;

    use super::{string_element_dofs_collection, ElementDofs};
    use crate::base::{Dof, Element, SampleParams};
    use gemlab::{mesh::Samples, shapes::GeoKind};

    #[test]
    fn new_handles_errors() {
        let p = SampleParams::param_rod();
        assert_eq!(
            ElementDofs::new(2, Element::Rod(p), GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        let p = SampleParams::param_beam();
        assert_eq!(
            ElementDofs::new(2, Element::Beam(p), GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        let p = SampleParams::param_solid();
        assert_eq!(
            ElementDofs::new(2, Element::Solid(p), GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        let p = SampleParams::param_porous_sld_liq();
        assert_eq!(
            ElementDofs::new(2, Element::PorousSldLiq(p), GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        let p = SampleParams::param_porous_sld_liq_gas();
        assert_eq!(
            ElementDofs::new(2, Element::PorousSldLiqGas(p), GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    fn new_works_2d() {
        let pa0 = SampleParams::param_diffusion();
        let pa = SampleParams::param_rod();
        let pb = SampleParams::param_beam();
        let pc = SampleParams::param_solid();
        let pd = SampleParams::param_porous_liq();
        let pe = SampleParams::param_porous_liq_gas();
        let pf = SampleParams::param_porous_sld_liq();
        let pg = SampleParams::param_porous_sld_liq_gas();
        let a0 = ElementDofs::new(2, Element::Diffusion(pa0), GeoKind::Tri3).unwrap();
        let a = ElementDofs::new(2, Element::Rod(pa), GeoKind::Lin2).unwrap();
        let b = ElementDofs::new(2, Element::Beam(pb), GeoKind::Lin2).unwrap();
        let c = ElementDofs::new(2, Element::Solid(pc), GeoKind::Tri3).unwrap();
        let d = ElementDofs::new(2, Element::PorousLiq(pd), GeoKind::Tri3).unwrap();
        let e = ElementDofs::new(2, Element::PorousLiqGas(pe), GeoKind::Tri3).unwrap();
        let f = ElementDofs::new(2, Element::PorousSldLiq(pf), GeoKind::Tri6).unwrap();
        let g = ElementDofs::new(2, Element::PorousSldLiqGas(pg), GeoKind::Tri6).unwrap();
        assert_eq!(a0.dof_equation_pairs, &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]);
        assert_eq!(
            a.dof_equation_pairs,
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(
            b.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Rz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Rz, 5)]
            ]
        );
        assert_eq!(
            c.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5)]
            ]
        );
        assert_eq!(d.dof_equation_pairs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            e.dof_equation_pairs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            f.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 12)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3), (Dof::Pl, 13)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5), (Dof::Pl, 14)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7)],
                vec![(Dof::Ux, 8), (Dof::Uy, 9)],
                vec![(Dof::Ux, 10), (Dof::Uy, 11)]
            ]
        );
        assert_eq!(
            g.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Pl, 12), (Dof::Pg, 15)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3), (Dof::Pl, 13), (Dof::Pg, 16)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5), (Dof::Pl, 14), (Dof::Pg, 17)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7)],
                vec![(Dof::Ux, 8), (Dof::Uy, 9)],
                vec![(Dof::Ux, 10), (Dof::Uy, 11)]
            ]
        );
    }

    #[test]
    fn new_works_3d() {
        let pa0 = SampleParams::param_diffusion();
        let pa = SampleParams::param_rod();
        let pb = SampleParams::param_beam();
        let pc = SampleParams::param_solid();
        let pd = SampleParams::param_porous_liq();
        let pe = SampleParams::param_porous_liq_gas();
        let pf = SampleParams::param_porous_sld_liq();
        let pg = SampleParams::param_porous_sld_liq_gas();
        let a0 = ElementDofs::new(3, Element::Diffusion(pa0), GeoKind::Tri3).unwrap();
        let a = ElementDofs::new(3, Element::Rod(pa), GeoKind::Lin2).unwrap();
        let b = ElementDofs::new(3, Element::Beam(pb), GeoKind::Lin2).unwrap();
        let c = ElementDofs::new(3, Element::Solid(pc), GeoKind::Tri3).unwrap();
        let d = ElementDofs::new(3, Element::PorousLiq(pd), GeoKind::Tri3).unwrap();
        let e = ElementDofs::new(3, Element::PorousLiqGas(pe), GeoKind::Tri3).unwrap();
        let f = ElementDofs::new(3, Element::PorousSldLiq(pf), GeoKind::Tri6).unwrap();
        let g = ElementDofs::new(3, Element::PorousSldLiqGas(pg), GeoKind::Tri6).unwrap();
        assert_eq!(a0.dof_equation_pairs, &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]);
        assert_eq!(
            a.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(
            b.dof_equation_pairs,
            vec![
                vec![
                    (Dof::Ux, 0),
                    (Dof::Uy, 1),
                    (Dof::Uz, 2),
                    (Dof::Rx, 3),
                    (Dof::Ry, 4),
                    (Dof::Rz, 5)
                ],
                vec![
                    (Dof::Ux, 6),
                    (Dof::Uy, 7),
                    (Dof::Uz, 8),
                    (Dof::Rx, 9),
                    (Dof::Ry, 10),
                    (Dof::Rz, 11)
                ]
            ]
        );
        assert_eq!(
            c.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8)]
            ]
        );
        assert_eq!(d.dof_equation_pairs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            e.dof_equation_pairs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            f.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2), (Dof::Pl, 18)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5), (Dof::Pl, 19)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8), (Dof::Pl, 20)],
                vec![(Dof::Ux, 9), (Dof::Uy, 10), (Dof::Uz, 11)],
                vec![(Dof::Ux, 12), (Dof::Uy, 13), (Dof::Uz, 14)],
                vec![(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Uz, 17)]
            ]
        );
        assert_eq!(
            g.dof_equation_pairs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2), (Dof::Pl, 18), (Dof::Pg, 21)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5), (Dof::Pl, 19), (Dof::Pg, 22)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8), (Dof::Pl, 20), (Dof::Pg, 23)],
                vec![(Dof::Ux, 9), (Dof::Uy, 10), (Dof::Uz, 11)],
                vec![(Dof::Ux, 12), (Dof::Uy, 13), (Dof::Uz, 14)],
                vec![(Dof::Ux, 15), (Dof::Uy, 16), (Dof::Uz, 17)]
            ]
        );
    }

    #[test]
    fn new_collection_handles_errors() {
        let mesh = Samples::one_tri6();
        let p2 = SampleParams::param_solid();
        let elements = HashMap::from([(2, Element::Solid(p2))]);
        assert_eq!(
            ElementDofs::new_collection(&mesh, &elements).err(),
            Some("cannot find CellAttributeId in elements map to create new ElementDofs collection")
        );
        let p1 = SampleParams::param_rod();
        let elements = HashMap::from([(1, Element::Rod(p1))]);
        assert_eq!(
            ElementDofs::new_collection(&mesh, &elements).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
    }

    #[test]
    fn display_works() {
        let p = SampleParams::param_porous_sld_liq();
        let ed = ElementDofs::new(1, Element::PorousSldLiq(p), GeoKind::Tri6).unwrap();
        assert_eq!(
            format!("{}", ed),
            "0: [(Ux, 0), (Uy, 1), (Pl, 12)]\n\
             1: [(Ux, 2), (Uy, 3), (Pl, 13)]\n\
             2: [(Ux, 4), (Uy, 5), (Pl, 14)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             5: [(Ux, 10), (Uy, 11)]\n\
             (Pl @ Some(12), Pg @ None, T @ None)\n"
        );
    }

    #[test]
    fn new_collection_works() {
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_beam();
        let elements = HashMap::from([
            (1, Element::PorousSldLiq(p1)),
            (2, Element::Solid(p2)),
            (3, Element::Beam(p3)),
        ]);
        let collection = ElementDofs::new_collection(&mesh, &elements).unwrap();
        let s = string_element_dofs_collection(&elements, &collection);
        assert_eq!(
            format!("{}", s),
            "Elements: DOFs and local equation numbers\n\
             =========================================\n\
             1 → PorousSldLiq → Qua8\n\
             0: [(Ux, 0), (Uy, 1), (Pl, 16)]\n\
             1: [(Ux, 2), (Uy, 3), (Pl, 17)]\n\
             2: [(Ux, 4), (Uy, 5), (Pl, 18)]\n\
             3: [(Ux, 6), (Uy, 7), (Pl, 19)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             5: [(Ux, 10), (Uy, 11)]\n\
             6: [(Ux, 12), (Uy, 13)]\n\
             7: [(Ux, 14), (Uy, 15)]\n\
             (Pl @ Some(16), Pg @ None, T @ None)\n\
             -----------------------------------------\n\
             2 → Solid → Tri6\n\
             0: [(Ux, 0), (Uy, 1)]\n\
             1: [(Ux, 2), (Uy, 3)]\n\
             2: [(Ux, 4), (Uy, 5)]\n\
             3: [(Ux, 6), (Uy, 7)]\n\
             4: [(Ux, 8), (Uy, 9)]\n\
             5: [(Ux, 10), (Uy, 11)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n\
             3 → Beam → Lin2\n\
             0: [(Ux, 0), (Uy, 1), (Rz, 2)]\n\
             1: [(Ux, 3), (Uy, 4), (Rz, 5)]\n\
             (Pl @ None, Pg @ None, T @ None)\n\
             -----------------------------------------\n"
        );
    }
}
