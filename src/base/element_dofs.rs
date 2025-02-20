use super::{Dof, Elem};
use crate::StrError;
use gemlab::mesh::GeoKind;
use serde::{Deserialize, Serialize};
use std::fmt;

/// Defines the allowed GeoKinds that can be used with PorousSld{...} elements
pub const POROUS_SLD_GEO_KIND_ALLOWED: [GeoKind; 7] = [
    // Tri
    GeoKind::Tri6,
    GeoKind::Tri15,
    // Qua
    GeoKind::Qua8,
    GeoKind::Qua9,
    GeoKind::Qua17,
    // Tet
    GeoKind::Tet10,
    // Hex
    GeoKind::Hex20,
];

/// Holds information of an (Element,GeoKind) pair such as DOFs and local equation numbers
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
#[derive(Clone, Debug, Deserialize, Serialize)]
pub struct ElementDofs {
    /// Holds all cell DOF keys and local equation numbers
    ///
    /// **Notes:** The outer array has length = nnode.
    /// The inner arrays have variable lengths = ndof at the node.
    /// The inner arrays contain pairs of Dof and local equation numbers.
    pub dofs: Vec<Vec<(Dof, usize)>>,

    /// Dimension of the local system of equations
    ///
    /// **Note:** This is equal to the total number of DOFs in the cell
    pub n_equation: usize,

    /// Local equation number of the first Dof::Pl
    pub eq_first_pl: Option<usize>,

    /// Local equation number of the first Dof::Pg
    pub eq_first_pg: Option<usize>,

    /// Local equation number of the first Dof::T
    pub eq_first_tt: Option<usize>,
}

impl ElementDofs {
    /// Allocates a new instance
    pub fn new(ndim: usize, element: Elem, kind: GeoKind) -> Result<Self, StrError> {
        // check
        let rod_or_beam = match element {
            Elem::Rod(..) => true,
            Elem::Beam(..) => true,
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
            Elem::Diffusion(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::T, count)); count += 1;
                }
            }
            Elem::Rod(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
            }
            Elem::Beam(..) => {
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
            Elem::Solid(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Ux, count)); count += 1;
                    dofs[m].push((Dof::Uy, count)); count += 1;
                    if ndim == 3 {
                        dofs[m].push((Dof::Uz, count)); count += 1;
                    }
                }
            }
            Elem::PorousLiq(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
            }
            Elem::PorousLiqGas(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                    dofs[m].push((Dof::Pg, count)); count += 1;
                }
            }
            Elem::PorousSldLiq(..) => {
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
            Elem::PorousSldLiqGas(..) => {
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
            dofs,
            n_equation: count,
            eq_first_pl,
            eq_first_pg,
            eq_first_tt,
        })
    }
}

impl fmt::Display for ElementDofs {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        for m in 0..self.dofs.len() {
            write!(f, "{}: {:?}\n", m, self.dofs[m]).unwrap();
        }
        write!(
            f,
            "(Pl @ {:?}, Pg @ {:?}, T @ {:?})\n",
            self.eq_first_pl, self.eq_first_pg, self.eq_first_tt
        )
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ElementDofs;
    use crate::base::{Dof, Elem, ParamDiffusion, ParamPorousLiq, ParamPorousLiqGas};
    use crate::base::{ParamBeam, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRod, ParamSolid};
    use gemlab::mesh::GeoKind;

    #[test]
    fn new_handles_errors() {
        let p = ParamRod::sample();
        assert_eq!(
            ElementDofs::new(2, Elem::Rod(p), GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        let p = ParamBeam::sample();
        assert_eq!(
            ElementDofs::new(2, Elem::Beam(p), GeoKind::Tri3).err(),
            Some("cannot set Rod or Beam with a non-Lin GeoClass")
        );
        let p = ParamSolid::sample_linear_elastic();
        assert_eq!(
            ElementDofs::new(2, Elem::Solid(p), GeoKind::Lin2).err(),
            Some("GeoClass::Lin is reserved for Rod or Beam")
        );
        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        assert_eq!(
            ElementDofs::new(2, Elem::PorousSldLiq(p), GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiq with given GeoKind")
        );
        let p = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        assert_eq!(
            ElementDofs::new(2, Elem::PorousSldLiqGas(p), GeoKind::Tri3).err(),
            Some("cannot set PorousSldLiqGas with given GeoKind")
        );
    }

    #[test]
    fn new_works_2d() {
        let pa = ParamDiffusion::sample();
        let pb = ParamRod::sample();
        let pc = ParamBeam::sample();
        let pd = ParamSolid::sample_linear_elastic();
        let pe = ParamPorousLiq::sample_brooks_corey_constant();
        let pf = ParamPorousLiqGas::sample_brooks_corey_constant();
        let pg = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let ph = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let a = ElementDofs::new(2, Elem::Diffusion(pa), GeoKind::Tri3).unwrap();
        let b = ElementDofs::new(2, Elem::Rod(pb), GeoKind::Lin2).unwrap();
        let c = ElementDofs::new(2, Elem::Beam(pc), GeoKind::Lin2).unwrap();
        let d = ElementDofs::new(2, Elem::Solid(pd), GeoKind::Tri3).unwrap();
        let e = ElementDofs::new(2, Elem::PorousLiq(pe), GeoKind::Tri3).unwrap();
        let f = ElementDofs::new(2, Elem::PorousLiqGas(pf), GeoKind::Tri3).unwrap();
        let g = ElementDofs::new(2, Elem::PorousSldLiq(pg), GeoKind::Tri6).unwrap();
        let h = ElementDofs::new(2, Elem::PorousSldLiqGas(ph), GeoKind::Tri6).unwrap();
        assert_eq!(a.dofs, &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]);
        assert_eq!(
            b.dofs,
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(
            c.dofs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Rz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Rz, 5)]
            ]
        );
        assert_eq!(
            d.dofs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1)],
                vec![(Dof::Ux, 2), (Dof::Uy, 3)],
                vec![(Dof::Ux, 4), (Dof::Uy, 5)]
            ]
        );
        assert_eq!(e.dofs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            f.dofs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            g.dofs,
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
            h.dofs,
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
        let pa = ParamDiffusion::sample();
        let pb = ParamRod::sample();
        let pc = ParamBeam::sample();
        let pd = ParamSolid::sample_linear_elastic();
        let pe = ParamPorousLiq::sample_brooks_corey_constant();
        let pf = ParamPorousLiqGas::sample_brooks_corey_constant();
        let pg = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let ph = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let a = ElementDofs::new(3, Elem::Diffusion(pa), GeoKind::Tri3).unwrap();
        let b = ElementDofs::new(3, Elem::Rod(pb), GeoKind::Lin2).unwrap();
        let c = ElementDofs::new(3, Elem::Beam(pc), GeoKind::Lin2).unwrap();
        let d = ElementDofs::new(3, Elem::Solid(pd), GeoKind::Tri3).unwrap();
        let e = ElementDofs::new(3, Elem::PorousLiq(pe), GeoKind::Tri3).unwrap();
        let f = ElementDofs::new(3, Elem::PorousLiqGas(pf), GeoKind::Tri3).unwrap();
        let g = ElementDofs::new(3, Elem::PorousSldLiq(pg), GeoKind::Tri6).unwrap();
        let h = ElementDofs::new(3, Elem::PorousSldLiqGas(ph), GeoKind::Tri6).unwrap();
        assert_eq!(a.dofs, &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]);
        assert_eq!(
            b.dofs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(
            c.dofs,
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
            d.dofs,
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)],
                vec![(Dof::Ux, 6), (Dof::Uy, 7), (Dof::Uz, 8)]
            ]
        );
        assert_eq!(e.dofs, &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]);
        assert_eq!(
            f.dofs,
            vec![
                vec![(Dof::Pl, 0), (Dof::Pg, 1)],
                vec![(Dof::Pl, 2), (Dof::Pg, 3)],
                vec![(Dof::Pl, 4), (Dof::Pg, 5)]
            ]
        );
        assert_eq!(
            g.dofs,
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
            h.dofs,
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
    fn display_works() {
        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let ed = ElementDofs::new(1, Elem::PorousSldLiq(p), GeoKind::Tri6).unwrap();
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
    fn derive_works() {
        let p1 = ParamSolid::sample_linear_elastic();
        let dofs = ElementDofs::new(2, Elem::Solid(p1), GeoKind::Tri3).unwrap();
        let clone = dofs.clone();
        let str_ori = format!("{:?}", dofs).to_string();
        assert_eq!(format!("{:?}", clone), str_ori);
        println!("{:?}", dofs);
        // serialize
        let json = serde_json::to_string(&dofs).unwrap();
        // deserialize
        let read: ElementDofs = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read), str_ori);
    }
}
