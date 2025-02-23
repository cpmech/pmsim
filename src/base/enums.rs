use super::{ParamBeam, ParamDiffusion, ParamRod};
use super::{ParamPorousLiq, ParamPorousLiqGas, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
use serde::{Deserialize, Serialize};

/// Defines degrees-of-freedom (DOF) types
///
/// Note: The fixed numbering scheme assists in sorting the DOFs.
#[derive(Clone, Copy, Debug, Eq, Hash, PartialEq, PartialOrd, Ord, Deserialize, Serialize)]
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

    /// Primary scalar quantity for diffusion problems (such as the temperature)
    ///
    /// Examples: Liquid pressure, temperature
    Phi = 6,

    /// Liquid pressure
    Pl = 7,

    /// Gas pressure
    Pg = 8,

    /// Free-surface-output (fso) enrichment
    Fso = 9,
}

/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Nbc {
    /// Normal distributed load
    Qn,

    /// Distributed load parallel to x
    Qx,

    /// Distributed load parallel to y
    Qy,

    /// Distributed load parallel to z
    Qz,

    /// Liquid flux
    Ql,

    /// Gas flux
    Qg,

    /// Heat flux
    Qt,

    /// Heat convection
    ///
    /// The value in parenthesis is constant and corresponds to is the convection coefficient `cc`.
    /// The specified value is the environment temperature `T∞`.
    Cv(f64),
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
            Nbc::Qn => solid(),
            Nbc::Qx => solid(),
            Nbc::Qy => solid(),
            Nbc::Qz => solid(),
            Nbc::Ql => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pl, count)); count += 1;
                }
            }
            Nbc::Qg => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Pg, count)); count += 1;
                }
            }
            Nbc::Qt => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Phi, count)); count += 1;
                }
            }
            Nbc::Cv(..) => {
                for m in 0..nnode {
                    dofs[m].push((Dof::Phi, count)); count += 1;
                }
            }
        }
        dofs
    }

    /// Indicates whether this NBC contributes to the Jacobian matrix or not
    pub fn contributes_to_jacobian_matrix(&self) -> bool {
        match self {
            Nbc::Qn => false,
            Nbc::Qx => false,
            Nbc::Qy => false,
            Nbc::Qz => false,
            Nbc::Ql => false,
            Nbc::Qg => false,
            Nbc::Qt => false,
            Nbc::Cv(..) => true,
        }
    }
}

/// Defines point boundary conditions (e.g., point loads)
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Pbc {
    /// Concentrated load parallel to x
    Fx,

    /// Concentrated load parallel to y
    Fy,

    /// Concentrated load parallel to z
    Fz,
}

impl Pbc {
    /// Returns the DOF corresponding to the concentrated load
    pub fn dof(&self) -> Dof {
        match self {
            Pbc::Fx => Dof::Ux,
            Pbc::Fy => Dof::Uy,
            Pbc::Fz => Dof::Uz,
        }
    }
}

/// Defines the kind of strain for geometrically non-linear analysis
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum GnlStrain {
    /// Engineering strain
    ///
    /// ```text
    ///      L - L₀
    /// ε  = ——————
    ///  ᴱ     L₀
    /// ```
    Eng,

    /// Green-Lagrange strain
    ///
    /// ```text
    ///      L² - L₀²
    /// ε  = ————————
    ///  ᴳ    2 L₀²
    /// ```
    Green,

    /// Logarithmic strain
    ///
    /// ```text
    ///         ⎛ L  ⎞
    /// ε  = ln ⎜————⎟
    ///  ᴸ      ⎝ L₀ ⎠
    /// ```
    Log,
}

/// Defines how stresses are initialized
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
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
#[derive(Clone, Copy, Debug, Deserialize, Serialize)]
pub enum Elem {
    Diffusion(ParamDiffusion),
    Rod(ParamRod),
    Beam(ParamBeam),
    Solid(ParamSolid),
    PorousLiq(ParamPorousLiq),
    PorousLiqGas(ParamPorousLiqGas),
    PorousSldLiq(ParamPorousSldLiq),
    PorousSldLiqGas(ParamPorousSldLiqGas),
}

impl Elem {
    /// Returns the name of the Element
    pub fn name(&self) -> String {
        match self {
            Elem::Diffusion(..) => "Diffusion".to_string(),
            Elem::Rod(..) => "Rod".to_string(),
            Elem::Beam(..) => "Beam".to_string(),
            Elem::Solid(..) => "Solid".to_string(),
            Elem::PorousLiq(..) => "PorousLiq".to_string(),
            Elem::PorousLiqGas(..) => "PorousLiqGas".to_string(),
            Elem::PorousSldLiq(..) => "PorousSldLiq".to_string(),
            Elem::PorousSldLiqGas(..) => "PorousSldLiqGas".to_string(),
        }
    }

    /// Returns the number of integration (Gauss) points
    pub fn ngauss(&self) -> Option<usize> {
        match self {
            Elem::Diffusion(param) => param.ngauss,
            Elem::Rod(param) => param.ngauss,
            Elem::Beam(param) => param.ngauss,
            Elem::Solid(param) => param.ngauss,
            Elem::PorousLiq(param) => param.ngauss,
            Elem::PorousLiqGas(param) => param.ngauss,
            Elem::PorousSldLiq(param) => param.ngauss,
            Elem::PorousSldLiqGas(param) => param.ngauss,
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{Dof, Elem, Init, Nbc, Pbc};
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

        // nbc
        let cv_ori = Nbc::Cv(1.2);
        let cv = cv_ori.clone();
        assert_eq!(format!("{:?}", cv), "Cv(1.2)");

        // pbc
        let fx_ori = Pbc::Fx;
        let fx = fx_ori.clone();
        assert_eq!(format!("{:?}", fx), "Fx");
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
        let e = Elem::Diffusion(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Diffusion");
        assert_eq!(e.ngauss(), None);

        let p = ParamRod::sample();
        let e = Elem::Rod(p);
        let e_clone = e.clone();
        assert_eq!(
            format!("{:?}", e),
            "Rod(ParamRod { gnl: None, density: 1.0, young: 1000.0, area: 1.0, ngauss: None })"
        );
        assert_eq!(format!("{}", e_clone.name()), "Rod");
        assert_eq!(e.ngauss(), None);

        let p = ParamBeam::sample();
        let e = Elem::Beam(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Beam");
        assert_eq!(e.ngauss(), None);

        let p = ParamSolid::sample_linear_elastic();
        let e = Elem::Solid(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "Solid");
        assert_eq!(e.ngauss(), None);

        let p = ParamPorousLiq::sample_brooks_corey_constant();
        let e = Elem::PorousLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiq");
        assert_eq!(e.ngauss(), None);

        let p = ParamPorousLiqGas::sample_brooks_corey_constant();
        let e = Elem::PorousLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousLiqGas");
        assert_eq!(e.ngauss(), None);

        let p = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let e = Elem::PorousSldLiq(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiq");
        assert_eq!(e.ngauss(), None);

        let p = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let e = Elem::PorousSldLiqGas(p);
        let e_clone = e.clone();
        assert_eq!(format!("{}", e_clone.name()), "PorousSldLiqGas");
        assert_eq!(e.ngauss(), None);
    }

    #[test]
    fn nbc_methods_work() {
        let qn = Nbc::Qn;
        assert_eq!(
            qn.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);

        let qx = Nbc::Qx;
        assert_eq!(
            qx.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);

        let qy = Nbc::Qy;
        assert_eq!(
            qy.dof_equation_pairs(2, 2),
            vec![vec![(Dof::Ux, 0), (Dof::Uy, 1)], vec![(Dof::Ux, 2), (Dof::Uy, 3)]]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);

        let qn = Nbc::Qn;
        assert_eq!(
            qn.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qn.contributes_to_jacobian_matrix(), false);

        let qx = Nbc::Qx;
        assert_eq!(
            qx.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qx.contributes_to_jacobian_matrix(), false);

        let qy = Nbc::Qy;
        assert_eq!(
            qy.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qy.contributes_to_jacobian_matrix(), false);

        let qz = Nbc::Qz;
        assert_eq!(
            qz.dof_equation_pairs(3, 2),
            vec![
                vec![(Dof::Ux, 0), (Dof::Uy, 1), (Dof::Uz, 2)],
                vec![(Dof::Ux, 3), (Dof::Uy, 4), (Dof::Uz, 5)]
            ]
        );
        assert_eq!(qz.contributes_to_jacobian_matrix(), false);

        let ql = Nbc::Ql;
        assert_eq!(
            ql.dof_equation_pairs(2, 3),
            &[[(Dof::Pl, 0)], [(Dof::Pl, 1)], [(Dof::Pl, 2)]]
        );
        assert_eq!(ql.contributes_to_jacobian_matrix(), false);

        let qt = Nbc::Qt;
        assert_eq!(
            qt.dof_equation_pairs(2, 3),
            &[[(Dof::Phi, 0)], [(Dof::Phi, 1)], [(Dof::Phi, 2)]]
        );
        assert_eq!(qt.contributes_to_jacobian_matrix(), false);

        let qg = Nbc::Qg;
        assert_eq!(
            qg.dof_equation_pairs(2, 3),
            &[[(Dof::Pg, 0)], [(Dof::Pg, 1)], [(Dof::Pg, 2)]]
        );
        assert_eq!(qg.contributes_to_jacobian_matrix(), false);

        let cv = Nbc::Cv(0.5);
        assert_eq!(
            cv.dof_equation_pairs(2, 3),
            &[[(Dof::Phi, 0)], [(Dof::Phi, 1)], [(Dof::Phi, 2)]]
        );
        assert_eq!(cv.contributes_to_jacobian_matrix(), true);
    }

    #[test]
    fn pbc_methods_work() {
        let fx = Pbc::Fx;
        assert_eq!(fx.dof(), Dof::Ux);

        let fy = Pbc::Fy;
        assert_eq!(fy.dof(), Dof::Uy);

        let fz = Pbc::Fz;
        assert_eq!(fz.dof(), Dof::Uz);
    }
}
