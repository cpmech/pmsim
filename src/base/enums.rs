use super::{
    ParamBeam, ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRod, ParamSolid,
};
use gemlab::shapes::GeoKind;
use std::fmt;

/// Defines a function to calculate boundary conditions values
///
/// This is a function of (t) where t is time
pub type FnBc = fn(t: f64) -> f64;

/// Defines degrees-of-freedom (DOF) types
///
/// Note: The fixed numbering scheme assists in sorting the DOFs.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord)]
pub enum Dof {
    /// Displacement along the first dimension
    Ux = 0,

    /// Displacement along the second dimension
    Uy = 1,

    /// Displacement along the third dimension
    Uz = 2,

    /// Rotation around the first axis
    Rx = 3,

    /// Rotation around the second axis
    Ry = 4,

    /// Rotation around the third axis
    Rz = 5,

    /// Temperature
    T = 6,

    /// Liquid pressure
    Pl = 7,

    /// Gas pressure
    Pg = 8,

    /// Free-surface-output (fso) enrichment
    Fso = 9,
}

/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy)]
pub enum Nbc {
    /// Normal distributed load
    Qn(FnBc),

    /// Distributed load parallel to x
    Qx(FnBc),

    /// Distributed load parallel to y
    Qy(FnBc),

    /// Distributed load parallel to z
    Qz(FnBc),

    /// Liquid flux
    Ql(FnBc),

    /// Gas flux
    Qg(FnBc),

    /// Convection
    ///
    /// The first value is the convection coefficient `cc`
    /// The second value is the environment temperature `T`
    Cv(f64, FnBc),
}

impl Nbc {
    /// Returns the boundary cell DOF keys and local equation numbers
    ///
    /// **Notes:** The outer array has length = nnode.
    /// The inner arrays have lengths = ndof at the node.
    #[rustfmt::skip]
    pub fn dof_equation_pairs(&self, ndim: usize, nnode: usize) -> Vec<Vec<(Dof, usize)>> {
        let mut dofs = vec![Vec::new(); nnode];
        let mut count = 0;
        let mut solid = || {
            for m in 0..nnode {
                dofs[m].push((Dof::Ux, count)); count += 1;
                dofs[m].push((Dof::Uy, count)); count += 1;
                if ndim == 3 {
                    dofs[m].push((Dof::Uz, count)); count += 1;
                }
            }
        };
        match self {
            Nbc::Qn(..) => solid(),
            Nbc::Qx(..) => solid(),
            Nbc::Qy(..) => solid(),
            Nbc::Qz(..) => solid(),
            Nbc::Ql(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
            }
            Nbc::Qg(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pg, count)); count += 1;
                }
            }
            Nbc::Cv(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::T, count)); count += 1;
                }
            }
        }
        dofs
    }

    /// Indicates whether this NBC contributes to the Jacobian matrix or not
    pub fn contributes_to_jacobian_matrix(&self) -> bool {
        match self {
            Nbc::Qn(..) => false,
            Nbc::Qx(..) => false,
            Nbc::Qy(..) => false,
            Nbc::Qz(..) => false,
            Nbc::Ql(..) => false,
            Nbc::Qg(..) => false,
            Nbc::Cv(..) => true,
        }
    }
}

impl fmt::Display for Nbc {
    fn fmt(&self, b: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Nbc::Qn(f) => write!(b, "Qn(0) = {:?}, Qn(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Qx(f) => write!(b, "Qx(0) = {:?}, Qx(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Qy(f) => write!(b, "Qy(0) = {:?}, Qy(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Qz(f) => write!(b, "Qz(0) = {:?}, Qz(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Ql(f) => write!(b, "Ql(0) = {:?}, Ql(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Qg(f) => write!(b, "Qg(0) = {:?}, Qg(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Cv(cc, f) => write!(b, "cc = {:?}, T(0) = {:?}, T(1) = {:?}", cc, f(0.0), f(1.0)).unwrap(),
        }
        Ok(())
    }
}

/// Defines point boundary conditions (e.g., point loads)
#[derive(Clone, Copy)]
pub enum Pbc {
    /// Concentrated load parallel to x
    Fx(FnBc),

    /// Concentrated load parallel to y
    Fy(FnBc),

    /// Concentrated load parallel to z
    Fz(FnBc),
}

impl Pbc {
    /// Returns the DOF corresponding to concentrated load
    pub fn dof(&self) -> Dof {
        match self {
            Pbc::Fx(..) => Dof::Ux,
            Pbc::Fy(..) => Dof::Uy,
            Pbc::Fz(..) => Dof::Uz,
        }
    }
}

impl fmt::Display for Pbc {
    fn fmt(&self, b: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Pbc::Fx(f) => write!(b, "Fx(0) = {:?}, Fx(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Pbc::Fy(f) => write!(b, "Fy(0) = {:?}, Fy(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Pbc::Fz(f) => write!(b, "Fz(0) = {:?}, Fz(1) = {:?}", f(0.0), f(1.0)).unwrap(),
        }
        Ok(())
    }
}

/// Defines how stresses are initialized
#[derive(Clone, Copy, Debug)]
pub enum Init {
    /// Geostatic initial state with data = (overburden,total_stress)
    ///
    /// # Note
    ///
    /// * The argument is the overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial isotropic stress state with σ_xx = σ_yy = σ_zz = value
    Isotropic(f64),

    /// Zero initial stress state
    Zero,
}

/// Defines the element type
#[derive(Clone, Copy, Debug)]
pub enum Element {
    Rod(ParamRod),
    Beam(ParamBeam),
    Solid(ParamSolid),
    PorousLiq(ParamPorousLiq),
    PorousLiqGas(ParamPorousLiqGas),
    PorousSldLiq(ParamPorousSldLiq),
    PorousSldLiqGas(ParamPorousSldLiqGas),
}

impl Element {
    /// Returns the name of the Element
    pub fn name(&self) -> String {
        match self {
            Element::Rod(..) => "Rod".to_string(),
            Element::Beam(..) => "Beam".to_string(),
            Element::Solid(..) => "Solid".to_string(),
            Element::PorousLiq(..) => "PorousLiq".to_string(),
            Element::PorousLiqGas(..) => "PorousLiqGas".to_string(),
            Element::PorousSldLiq(..) => "PorousSldLiq".to_string(),
            Element::PorousSldLiqGas(..) => "PorousSldLiqGas".to_string(),
        }
    }
}

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

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use crate::base::SampleParams;

    use super::{Dof, Element, Init, Nbc, Pbc};
    use std::{cmp::Ordering, collections::HashSet};

    #[test]
    fn dof_derive_works() {
        let ux = Dof::Ux;
        let ux_clone = ux.clone();
        assert_eq!(format!("{:?}", ux), "Ux");
        assert_eq!(ux, ux_clone);

        let uy = Dof::Uy;
        assert!(ux < uy);
        assert_eq!(ux.cmp(&uy), Ordering::Less);

        let mut set = HashSet::new();
        set.insert(ux);
        assert_eq!(set.len(), 1);
    }

    #[test]
    fn nbc_and_pbc_deriv_works() {
        let qn_ori = Nbc::Qn(|t| -10.0 * (1.0 + t));
        let qn = qn_ori.clone();
        assert_eq!(format!("{}", qn), "Qn(0) = -10.0, Qn(1) = -20.0");

        let qx_ori = Nbc::Qx(|t| -20.0 * (1.0 + t));
        let qx = qx_ori;
        assert_eq!(format!("{}", qx), "Qx(0) = -20.0, Qx(1) = -40.0");

        let qy = Nbc::Qy(|t| -20.0 * (1.0 + t));
        assert_eq!(format!("{}", qy), "Qy(0) = -20.0, Qy(1) = -40.0");

        let qz = Nbc::Qz(|t| -20.0 * (1.0 + t));
        assert_eq!(format!("{}", qz), "Qz(0) = -20.0, Qz(1) = -40.0");

        let ql = Nbc::Ql(|t| -20.0 * (1.0 + t));
        assert_eq!(format!("{}", ql), "Ql(0) = -20.0, Ql(1) = -40.0");

        let qg = Nbc::Qg(|t| -20.0 * (1.0 + t));
        assert_eq!(format!("{}", qg), "Qg(0) = -20.0, Qg(1) = -40.0");

        let cv = Nbc::Cv(0.5, |t| 50.0 * (1.0 + t));
        assert_eq!(format!("{}", cv), "cc = 0.5, T(0) = 50.0, T(1) = 100.0");

        let fx_ori = Pbc::Fx(|t| -1.0 * (1.0 + t));
        let fx = fx_ori.clone();
        assert_eq!(format!("{}", fx), "Fx(0) = -1.0, Fx(1) = -2.0");

        let fy_ori = Pbc::Fy(|t| -1.0 * (1.0 + t));
        let fy = fy_ori;
        assert_eq!(format!("{}", fy), "Fy(0) = -1.0, Fy(1) = -2.0");

        let fz = Pbc::Fz(|t| -1.0 * (1.0 + t));
        assert_eq!(format!("{}", fz), "Fz(0) = -1.0, Fz(1) = -2.0");
    }

    #[test]
    fn init_derive_works() {
        let init = Init::Geostatic(123.456);
        let init_clone = init.clone();
        assert_eq!(format!("{:?}", init_clone), format!("{:?}", init));
    }

    #[test]
    fn element_derive_works() {
        let p = SampleParams::param_rod();
        let e = Element::Rod(p);
        let e_clone = e.clone();
        assert_eq!(
            format!("{:?}", e),
            "Rod(ParamRod { density: 2.0, young: 1000.0, area: 1.0 })"
        );
        assert_eq!(format!("{}", e_clone.name()), "Rod");

        let p = SampleParams::param_beam();
        let e = Element::Beam(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Beam");

        let p = SampleParams::param_solid();
        let e = Element::Solid(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Solid");

        let p = SampleParams::param_porous_liq();
        let e = Element::PorousLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiq");

        let p = SampleParams::param_porous_liq_gas();
        let e = Element::PorousLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiqGas");

        let p = SampleParams::param_porous_sld_liq();
        let e = Element::PorousSldLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiq");

        let p = SampleParams::param_porous_sld_liq_gas();
        let e = Element::PorousSldLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiqGas");
    }

    #[test]
    fn nbc_methods_work() {
        let qn = Nbc::Qn(|_| -10.0);
        assert_eq!(
            qn.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qn), "Qn(0) = -10.0, Qn(1) = -10.0");

        let qx = Nbc::Qx(|_| -10.0);
        assert_eq!(
            qx.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qx), "Qx(0) = -10.0, Qx(1) = -10.0");

        let qy = Nbc::Qy(|_| -10.0);
        assert_eq!(
            qy.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qy), "Qy(0) = -10.0, Qy(1) = -10.0");

        let qn = Nbc::Qn(|_| -10.0);
        assert_eq!(
            qn.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qn), "Qn(0) = -10.0, Qn(1) = -10.0");

        let qx = Nbc::Qx(|_| -10.0);
        assert_eq!(
            qx.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qx), "Qx(0) = -10.0, Qx(1) = -10.0");

        let qy = Nbc::Qy(|_| -10.0);
        assert_eq!(
            qy.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qy), "Qy(0) = -10.0, Qy(1) = -10.0");

        let qz = Nbc::Qz(|_| -10.0);
        assert_eq!(
            qz.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qz.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qz), "Qz(0) = -10.0, Qz(1) = -10.0");

        let ql = Nbc::Ql(|_| -10.0);
        assert_eq!(
            ql.dof_equation_pairs(2, 3),
            &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]
        );
        assert_eq!(ql.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", ql), "Ql(0) = -10.0, Ql(1) = -10.0");

        let qg = Nbc::Qg(|_| -10.0);
        assert_eq!(
            qg.dof_equation_pairs(2, 3),
            &[[(Dof::Pg, 0)], [(Dof::Pg, 1)], [(Dof::Pg, 2)]]
        );
        assert_eq!(qg.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qg), "Qg(0) = -10.0, Qg(1) = -10.0");

        let cv = Nbc::Cv(0.5, |_| 25.0);
        assert_eq!(
            cv.dof_equation_pairs(2, 3),
            &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]
        );
        assert_eq!(cv.contributes_to_jacobian_matrix(), true);
        assert_eq!(format!("{}", cv), "cc = 0.5, T(0) = 25.0, T(1) = 25.0");
    }

    #[test]
    fn pbc_dof_works() {
        let fx = Pbc::Fx(|_| -1.0);
        assert_eq!(fx.dof(), Dof::Ux);
        assert_eq!(format!("{}", fx), "Fx(0) = -1.0, Fx(1) = -1.0");

        let fy = Pbc::Fy(|_| -1.0);
        assert_eq!(fy.dof(), Dof::Uy);
        assert_eq!(format!("{}", fy), "Fy(0) = -1.0, Fy(1) = -1.0");

        let fz = Pbc::Fz(|_| -1.0);
        assert_eq!(fz.dof(), Dof::Uz);
        assert_eq!(format!("{}", fz), "Fz(0) = -1.0, Fz(1) = -1.0");
    }
}
