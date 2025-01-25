use super::{ParamBeam, ParamDiffusion, ParamRod};
use super::{ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
use serde::{Deserialize, Serialize};
use std::fmt;

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
    Ux(f64),

    /// Displacement along the second dimension
    Uy(f64),

    /// Displacement along the third dimension
    Uz(f64),

    /// Rotation around the first axis
    Rx(f64),

    /// Rotation around the second axis
    Ry(f64),

    /// Rotation around the third axis
    Rz(f64),

    /// Temperature
    T(f64),

    /// Liquid pressure
    Pl(f64),

    /// Gas pressure
    Pg(f64),
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
            Ebc::Ux(v) => write!(b, "Ux = {:?}", v).unwrap(),
            Ebc::Uy(v) => write!(b, "Uy = {:?}", v).unwrap(),
            Ebc::Uz(v) => write!(b, "Uz = {:?}", v).unwrap(),
            Ebc::Rx(v) => write!(b, "Rx = {:?}", v).unwrap(),
            Ebc::Ry(v) => write!(b, "Ry = {:?}", v).unwrap(),
            Ebc::Rz(v) => write!(b, "Rz = {:?}", v).unwrap(),
            Ebc::T(v) => write!(b, "T = {:?}", v).unwrap(),
            Ebc::Pl(v) => write!(b, "Pl = {:?}", v).unwrap(),
            Ebc::Pg(v) => write!(b, "Pg = {:?}", v).unwrap(),
        }
        Ok(())
    }
}

/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy)]
pub enum Nbc {
    /// Normal distributed load
    Qn(f64),

    /// Distributed load parallel to x
    Qx(f64),

    /// Distributed load parallel to y
    Qy(f64),

    /// Distributed load parallel to z
    Qz(f64),

    /// Liquid flux
    Ql(f64),

    /// Gas flux
    Qg(f64),

    /// Temperature flux
    Qt(f64),

    /// Temperature convection
    ///
    /// The first value is the convection coefficient `cc`
    /// The second value is the environment temperature `T_env`
    Cv(f64, f64),
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
            Nbc::Qn(v) => write!(b, "Qn = {:?}", v).unwrap(),
            Nbc::Qx(v) => write!(b, "Qx = {:?}", v).unwrap(),
            Nbc::Qy(v) => write!(b, "Qy = {:?}", v).unwrap(),
            Nbc::Qz(v) => write!(b, "Qz = {:?}", v).unwrap(),
            Nbc::Ql(v) => write!(b, "Ql = {:?}", v).unwrap(),
            Nbc::Qg(v) => write!(b, "Qg = {:?}", v).unwrap(),
            Nbc::Qt(v) => write!(b, "Qt = {:?}", v).unwrap(),
            Nbc::Cv(cc, tt_env) => write!(b, "cc = {:?}, T_env = {:?}", cc, tt_env).unwrap(),
        }
        Ok(())
    }
}

/// Defines point boundary conditions (e.g., point loads)
#[derive(Clone, Copy)]
pub enum Pbc {
    /// Concentrated load parallel to x
    Fx(f64),

    /// Concentrated load parallel to y
    Fy(f64),

    /// Concentrated load parallel to z
    Fz(f64),
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
            Pbc::Fx(v) => write!(b, "Fx = {:?}", v).unwrap(),
            Pbc::Fy(v) => write!(b, "Fy = {:?}", v).unwrap(),
            Pbc::Fz(v) => write!(b, "Fz = {:?}", v).unwrap(),
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
        let ebc_ux_ori = Ebc::Ux(10.0);
        let ebc_ux = ebc_ux_ori.clone();
        assert_eq!(format!("{}", ebc_ux), "Ux = 10.0");

        // nbc
        let qn_ori = Nbc::Qn(-10.0);
        let qn = qn_ori.clone();
        assert_eq!(format!("{}", qn), "Qn = -10.0");

        // pbc
        let fx_ori = Pbc::Fx(-1.0);
        let fx = fx_ori.clone();
        assert_eq!(format!("{}", fx), "Fx = -1.0");
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
        let ux = Ebc::Ux(20.0);
        assert_eq!(ux.dof(), Dof::Ux);
        assert_eq!(format!("{}", ux), "Ux = 20.0");

        let uy = Ebc::Uy(20.0);
        assert_eq!(uy.dof(), Dof::Uy);
        assert_eq!(format!("{}", uy), "Uy = 20.0");

        let uz = Ebc::Uz(20.0);
        assert_eq!(uz.dof(), Dof::Uz);
        assert_eq!(format!("{}", uz), "Uz = 20.0");

        let rx = Ebc::Rx(20.0);
        assert_eq!(rx.dof(), Dof::Rx);
        assert_eq!(format!("{}", rx), "Rx = 20.0");

        let ry = Ebc::Ry(20.0);
        assert_eq!(ry.dof(), Dof::Ry);
        assert_eq!(format!("{}", ry), "Ry = 20.0");

        let rz = Ebc::Rz(20.0);
        assert_eq!(rz.dof(), Dof::Rz);
        assert_eq!(format!("{}", rz), "Rz = 20.0");

        let pl = Ebc::Pl(20.0);
        assert_eq!(pl.dof(), Dof::Pl);
        assert_eq!(format!("{}", pl), "Pl = 20.0");

        let pg = Ebc::Pg(20.0);
        assert_eq!(pg.dof(), Dof::Pg);
        assert_eq!(format!("{}", pg), "Pg = 20.0");

        let tt = Ebc::T(20.0);
        assert_eq!(tt.dof(), Dof::T);
        assert_eq!(format!("{}", tt), "T = 20.0");
    }

    #[test]
    fn nbc_methods_work() {
        let qn = Nbc::Qn(-10.0);
        assert_eq!(
            qn.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qn), "Qn = -10.0");

        let qx = Nbc::Qx(-10.0);
        assert_eq!(
            qx.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qx), "Qx = -10.0");

        let qy = Nbc::Qy(-10.0);
        assert_eq!(
            qy.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qy), "Qy = -10.0");

        let qn = Nbc::Qn(-10.0);
        assert_eq!(
            qn.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qn), "Qn = -10.0");

        let qx = Nbc::Qx(-10.0);
        assert_eq!(
            qx.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qx), "Qx = -10.0");

        let qy = Nbc::Qy(-10.0);
        assert_eq!(
            qy.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qy), "Qy = -10.0");

        let qz = Nbc::Qz(-10.0);
        assert_eq!(
            qz.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qz.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qz), "Qz = -10.0");

        let ql = Nbc::Ql(-10.0);
        assert_eq!(
            ql.dof_equation_pairs(2, 3),
            &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]
        );
        assert_eq!(ql.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", ql), "Ql = -10.0");

        let qt = Nbc::Qt(-10.0);
        assert_eq!(
            qt.dof_equation_pairs(2, 3),
            &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]
        );
        assert_eq!(qt.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qt), "Qt = -10.0");

        let qg = Nbc::Qg(-10.0);
        assert_eq!(
            qg.dof_equation_pairs(2, 3),
            &[[(Dof::Pg, 0)], [(Dof::Pg, 1)], [(Dof::Pg, 2)]]
        );
        assert_eq!(qg.contributes_to_jacobian_matrix(), false);
        assert_eq!(format!("{}", qg), "Qg = -10.0");

        let cv = Nbc::Cv(0.5, 25.0);
        assert_eq!(
            cv.dof_equation_pairs(2, 3),
            &[[(Dof::T, 0)], [(Dof::T, 1)], [(Dof::T, 2)]]
        );
        assert_eq!(cv.contributes_to_jacobian_matrix(), true);
        assert_eq!(format!("{}", cv), "cc = 0.5, T_env = 25.0");
    }

    #[test]
    fn pbc_methods_work() {
        let fx = Pbc::Fx(-1.0);
        assert_eq!(fx.dof(), Dof::Ux);
        assert_eq!(format!("{}", fx), "Fx = -1.0");

        let fy = Pbc::Fy(-1.0);
        assert_eq!(fy.dof(), Dof::Uy);
        assert_eq!(format!("{}", fy), "Fy = -1.0");

        let fz = Pbc::Fz(-1.0);
        assert_eq!(fz.dof(), Dof::Uz);
        assert_eq!(format!("{}", fz), "Fz = -1.0");
    }
}
