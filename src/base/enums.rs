use super::{ParamBeam, ParamDiffusion, ParamRod};
use super::{ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
use serde::{Deserialize, Serialize};
use std::fmt;

/// Defines a function of time that returns f64 (e.g., to calculate boundary condition values)
pub type FnTime = fn(t: f64) -> f64;

/// Defines degrees-of-freedom (DOF) types
///
/// Note: The fixed numbering scheme assists in sorting the DOFs.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord, Serialize, Deserialize)]
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

/// Defines essential boundary conditions (EBC)
#[derive(Clone, Copy)]
pub enum Ebc {
    /// Displacement along the first dimension
    Ux(FnTime),

    /// Displacement along the second dimension
    Uy(FnTime),

    /// Displacement along the third dimension
    Uz(FnTime),

    /// Rotation around the first axis
    Rx(FnTime),

    /// Rotation around the second axis
    Ry(FnTime),

    /// Rotation around the third axis
    Rz(FnTime),

    /// Temperature
    T(FnTime),

    /// Liquid pressure
    Pl(FnTime),

    /// Gas pressure
    Pg(FnTime),
}

impl Ebc {
    /// Returns the DOF corresponding to the essential boundary condition
    pub fn dof(&self) -> Dof {
        match self {
            Ebc::Ux(..) => Dof::Ux,
            Ebc::Uy(..) => Dof::Uy,
            Ebc::Uz(..) => Dof::Uz,
            Ebc::Rx(..) => Dof::Rx,
            Ebc::Ry(..) => Dof::Ry,
            Ebc::Rz(..) => Dof::Rz,
            Ebc::T(..) => Dof::T,
            Ebc::Pl(..) => Dof::Pl,
            Ebc::Pg(..) => Dof::Pg,
        }
    }
}

impl fmt::Display for Ebc {
    fn fmt(&self, b: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Ebc::Ux(f) => write!(b, "Ux(0) = {:?}, Ux(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Uy(f) => write!(b, "Uy(0) = {:?}, Uy(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Uz(f) => write!(b, "Uz(0) = {:?}, Uz(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Rx(f) => write!(b, "Rx(0) = {:?}, Rx(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Ry(f) => write!(b, "Ry(0) = {:?}, Ry(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Rz(f) => write!(b, "Rz(0) = {:?}, Rz(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::T(f) => write!(b, "T(0) = {:?}, T(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Pl(f) => write!(b, "Pl(0) = {:?}, Pl(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Ebc::Pg(f) => write!(b, "Pg(0) = {:?}, Pg(1) = {:?}", f(0.0), f(1.0)).unwrap(),
        }
        Ok(())
    }
}

/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy)]
pub enum Nbc {
    /// Normal distributed load
    Qn(FnTime),

    /// Distributed load parallel to x
    Qx(FnTime),

    /// Distributed load parallel to y
    Qy(FnTime),

    /// Distributed load parallel to z
    Qz(FnTime),

    /// Liquid flux
    Ql(FnTime),

    /// Gas flux
    Qg(FnTime),

    /// Temperature flux
    Qt(FnTime),

    /// Temperature convection
    ///
    /// The first value is the convection coefficient `cc`
    /// The second value is the environment temperature `T`
    Cv(f64, FnTime),
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
            Nbc::Qt(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::T, count)); count += 1;
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
            Nbc::Qt(..) => false,
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
            Nbc::Qt(f) => write!(b, "Qt(0) = {:?}, Qt(1) = {:?}", f(0.0), f(1.0)).unwrap(),
            Nbc::Cv(cc, f) => write!(b, "cc = {:?}, T(0) = {:?}, T(1) = {:?}", cc, f(0.0), f(1.0)).unwrap(),
        }
        Ok(())
    }
}

/// Defines point boundary conditions (e.g., point loads)
#[derive(Clone, Copy)]
pub enum Pbc {
    /// Concentrated load parallel to x
    Fx(FnTime),

    /// Concentrated load parallel to y
    Fy(FnTime),

    /// Concentrated load parallel to z
    Fz(FnTime),
}

impl Pbc {
    /// Returns the DOF corresponding to the concentrated load
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
    /// Geostatic initial state with data = (overburden)
    ///
    /// # Note
    ///
    /// * The argument is the overburden stress (negative means compression) at the whole surface (z=z_max=height)
    /// * The datum is at y=0.0 (2D) or z=0.0 (3D)
    /// * The water table is at y=y_max=height (2D) or z=z_max=height (3D), thus only fully water-saturated states are considered
    Geostatic(f64),

    /// Initial homogeneous AND isotropic stress state with σ_xx = σ_yy = σ_zz = value
    ///
    /// All points in the mesh will be set with the same isotropic stress
    Isotropic(f64),

    /// Zero initial stress state
    Zero,
}

/// Defines the element type
#[derive(Clone, Copy, Debug)]
pub enum Etype {
    Diffusion(ParamDiffusion),
    Rod(ParamRod),
    Beam(ParamBeam),
    Solid(ParamSolid),
    PorousLiq(ParamPorousLiq),
    PorousLiqGas(ParamPorousLiqGas),
    PorousSldLiq(ParamPorousSldLiq),
    PorousSldLiqGas(ParamPorousSldLiqGas),
}

impl Etype {
    /// Returns the name of the Element
    pub fn name(&self) -> String {
        match self {
            Etype::Diffusion(..) => "Diffusion".to_string(),
            Etype::Rod(..) => "Rod".to_string(),
            Etype::Beam(..) => "Beam".to_string(),
            Etype::Solid(..) => "Solid".to_string(),
            Etype::PorousLiq(..) => "PorousLiq".to_string(),
            Etype::PorousLiqGas(..) => "PorousLiqGas".to_string(),
            Etype::PorousSldLiq(..) => "PorousSldLiq".to_string(),
            Etype::PorousSldLiqGas(..) => "PorousSldLiqGas".to_string(),
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Dof, Ebc, Etype, Init, Nbc, Pbc};
    use crate::base::{ParamBeam, ParamDiffusion, ParamPorousLiq, ParamPorousLiqGas};
    use crate::base::{ParamPorousSldLiq, ParamPorousSldLiqGas, ParamRod, ParamSolid};
    use std::{cmp::Ordering, collections::HashSet};

    #[test]
    fn dof_ebc_nbc_pbc_derives_work() {
        // dof
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

        // ebc
        let ebc_ux_ori = Ebc::Ux(|_| 10.0);
        let ebc_ux = ebc_ux_ori.clone();
        assert_eq!(format!("{}", ebc_ux), "Ux(0) = 10.0, Ux(1) = 10.0");

        // nbc
        let qn_ori = Nbc::Qn(|t| -10.0 * (1.0 + t));
        let qn = qn_ori.clone();
        assert_eq!(format!("{}", qn), "Qn(0) = -10.0, Qn(1) = -20.0");

        // pbc
        let fx_ori = Pbc::Fx(|t| -1.0 * (1.0 + t));
        let fx = fx_ori.clone();
        assert_eq!(format!("{}", fx), "Fx(0) = -1.0, Fx(1) = -2.0");
    }

    #[test]
    fn init_derive_works() {
        let init = Init::Geostatic(123.456);
        let init_clone = init.clone();
        assert_eq!(format!("{:?}", init_clone), format!("{:?}", init));
    }

    #[test]
    fn element_derive_works() {
        let p = ParamDiffusion::sample();
        let e = Etype::Diffusion(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Diffusion");

        let p = ParamRod::sample();
        let e = Etype::Rod(p);
        let e_clone = e.clone();
        assert_eq!(
            format!("{:?}", e),
            "Rod(ParamRod { density: 1.0, young: 1000.0, area: 1.0 })"
        );
        assert_eq!(format!("{}", e_clone.name()), "Rod");

        let p = ParamBeam::sample();
        let e = Etype::Beam(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Beam");

        let p = ParamSolid::sample_linear_elastic();
        let e = Etype::Solid(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Solid");

        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let e = Etype::PorousLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiq");

        let p = ParamPorousLiqGas::sample_brooks_corey_constant();
        let e = Etype::PorousLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiqGas");

        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let e = Etype::PorousSldLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiq");

        let p = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let e = Etype::PorousSldLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiqGas");
    }

    #[test]
    fn ebc_methods_work() {
        let ux = Ebc::Ux(|t| 10.0 + 10.0 * t);
        assert_eq!(ux.dof(), Dof::Ux);
        assert_eq!(format!("{}", ux), "Ux(0) = 10.0, Ux(1) = 20.0");

        let uy = Ebc::Uy(|t| 10.0 + 10.0 * t);
        assert_eq!(uy.dof(), Dof::Uy);
        assert_eq!(format!("{}", uy), "Uy(0) = 10.0, Uy(1) = 20.0");

        let uz = Ebc::Uz(|t| 10.0 + 10.0 * t);
        assert_eq!(uz.dof(), Dof::Uz);
        assert_eq!(format!("{}", uz), "Uz(0) = 10.0, Uz(1) = 20.0");

        let rx = Ebc::Rx(|t| 10.0 + 10.0 * t);
        assert_eq!(rx.dof(), Dof::Rx);
        assert_eq!(format!("{}", rx), "Rx(0) = 10.0, Rx(1) = 20.0");

        let ry = Ebc::Ry(|t| 10.0 + 10.0 * t);
        assert_eq!(ry.dof(), Dof::Ry);
        assert_eq!(format!("{}", ry), "Ry(0) = 10.0, Ry(1) = 20.0");

        let rz = Ebc::Rz(|t| 10.0 + 10.0 * t);
        assert_eq!(rz.dof(), Dof::Rz);
        assert_eq!(format!("{}", rz), "Rz(0) = 10.0, Rz(1) = 20.0");

        let pl = Ebc::Pl(|t| 10.0 + 10.0 * t);
        assert_eq!(pl.dof(), Dof::Pl);
        assert_eq!(format!("{}", pl), "Pl(0) = 10.0, Pl(1) = 20.0");

        let pg = Ebc::Pg(|t| 10.0 + 10.0 * t);
        assert_eq!(pg.dof(), Dof::Pg);
        assert_eq!(format!("{}", pg), "Pg(0) = 10.0, Pg(1) = 20.0");

        let tt = Ebc::T(|t| 10.0 + 10.0 * t);
        assert_eq!(tt.dof(), Dof::T);
        assert_eq!(format!("{}", tt), "T(0) = 10.0, T(1) = 20.0");
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

        let qt = Nbc::Qt(|_| -10.0);
        assert_eq!(
            qt.dof_equation_pairs(2, 3),
            &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]
        );
        assert_eq!(qt.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qt), "Qt(0) = -10.0, Qt(1) = -10.0");

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
    fn pbc_methods_work() {
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
