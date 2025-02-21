use super::{FemBase, FemState};
use crate::base::Essential;
use crate::StrError;
use russell_lab::Vector;
use russell_sparse::CooMatrix;

/// Assists in calculating essential (Dirichlet) boundary conditions (aka prescribed BCs)
pub struct BcPrescribed<'a> {
    /// All constant values
    ///
    /// (n_prescribed)
    constants: Vec<f64>,

    /// All multiplier functions
    ///
    /// (n_prescribed)
    multipliers: Vec<Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>>,

    /// Array with only the numbers of the prescribed DOFs
    ///
    /// (n_prescribed)
    pub equations: Vec<usize>,

    /// An array indicating which DOFs are prescribed
    ///
    /// (ndof; the total number of DOFs)
    pub flags: Vec<bool>,
}

impl<'a> BcPrescribed<'a> {
    /// Allocates a new instance
    pub fn new(base: &FemBase, essential: &'a Essential) -> Result<Self, StrError> {
        let n_prescribed = essential.size();
        let mut constants = Vec::with_capacity(n_prescribed);
        let mut multipliers = Vec::with_capacity(n_prescribed);
        let mut equations = Vec::with_capacity(n_prescribed);
        let mut flags = vec![false; base.dofs.size()];
        for (point_id, dof) in essential.keys() {
            let eq = base.dofs.eq(*point_id, *dof)?;
            let (constant, multiplier) = essential.get(*point_id, *dof);
            constants.push(constant);
            multipliers.push(multiplier);
            flags[eq] = true;
            equations.push(eq);
        }
        Ok(BcPrescribed {
            constants,
            multipliers,
            flags,
            equations,
        })
    }

    /// Returns the value of the prescribed DOF at given time
    ///
    /// The BC value is computed as follows:
    ///
    /// ```text
    /// value = constant * multiplier(t)
    /// ```
    pub fn value(&self, eq: usize, time: f64) -> f64 {
        match self.multipliers[eq] {
            Some(m) => self.constants[eq] * (m)(time),
            None => self.constants[eq],
        }
    }

    /// Returns the number of prescribed values
    pub fn size(&self) -> usize {
        self.constants.len()
    }

    /// Returns true if there is at least one non-zero constant
    pub fn has_non_zero(&self) -> bool {
        for constant in &self.constants {
            if *constant != 0.0 {
                return true;
            }
        }
        false
    }

    /// Assembles the contribution due to the prescribed DOFs into the global R vector (LMM)
    ///
    /// **LMM** means Lagrange Multiplier Method
    ///
    /// This function adds `Aᵀλ` to the global R vector at the non-prescribed equations and
    /// **sets** the prescribed equations to `A u - c`. Here, `c` is the prescribed value.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌         ┐
    ///  │  K   Aᵀ │ │ -δu │   │ R + Aᵀλ │
    ///  │         │ │     │ = │         │
    ///  │  A   0  │ │ -δλ │   │ A u - c │
    ///  └         ┘ └     ┘   └         ┘
    /// ```
    pub fn assemble_rr_lmm(&self, rr: &mut Vector, state: &FemState) {
        let ndof = self.flags.len();
        for p in 0..self.equations.len() {
            let i = self.equations[p];
            let j = ndof + p;
            let lambda = state.u[j];
            let c = self.value(p, state.t);
            rr[i] += lambda; // Aᵀ λ  →  1 * λ
            rr[j] = state.u[i] - c; // A u - c  →  1 * u - c
        }
    }

    /// Assembles the constraint matrix into the global K matrix (LMM)
    ///
    /// **LMM** means Lagrange Multiplier Method
    ///
    /// This function adds the constraints matrix (Aᵀ and A) to K.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌         ┐
    ///  │  K   Aᵀ │ │ -δu │   │ R + Aᵀλ │
    ///  │         │ │     │ = │         │
    ///  │  A   0  │ │ -δλ │   │ A u - c │
    ///  └         ┘ └     ┘   └         ┘
    /// ```
    pub fn assemble_kk_lmm(&self, kk: &mut CooMatrix) {
        let ndof = self.flags.len();
        for p in 0..self.equations.len() {
            let i = self.equations[p];
            let j = ndof + p;
            kk.put(i, j, 1.0).unwrap(); // Aᵀ
            kk.put(j, i, 1.0).unwrap(); // A
        }
    }

    /// Updates the diagonal of the global K matrix (RSM)
    ///
    /// **RSM** means Reduced-System Method
    ///
    /// This function put ones on the diagonal entries corresponding to the prescribed DOFs.
    ///
    /// The global system is symbolized by:
    ///
    /// ```text
    ///  ┌         ┐ ┌     ┐   ┌   ┐
    ///  │  K   0  │ │ -δu │   │ R │
    ///  │         │ │     │ = │   │
    ///  │  0   1  │ │  0  │   │ 0 │
    ///  └         ┘ └     ┘   └   ┘
    /// ```
    ///
    /// Note that the prescribed values are zero (homogeneous BCs).
    pub fn assemble_kk_rsm(&self, kk: &mut CooMatrix) {
        for eq in &self.equations {
            kk.put(*eq, *eq, 1.0).unwrap();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcPrescribed;
    use crate::base::{Dof, Elem, Essential, ParamBeam, ParamDiffusion};
    use crate::base::{ParamPorousLiq, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
    use crate::fem::FemBase;
    use gemlab::mesh::{Cell, GeoKind, Mesh, Point, Samples};

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

        let mut essential = Essential::new();
        essential.points(&[100], Dof::Ux, 0.0);
        assert_eq!(
            BcPrescribed::new(&base, &essential).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );

        let mut essential = Essential::new();
        essential.points(&[0], Dof::T, 0.0);
        assert_eq!(
            BcPrescribed::new(&base, &essential).err(),
            Some("cannot find the number of a (PointId, DOF) pair")
        );
    }

    #[test]
    fn bc_prescribed_array_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let mut essential = Essential::new();
        essential.points(&[0], Dof::T, 110.0);
        let array = BcPrescribed::new(&base, &essential).unwrap();
        assert_eq!(array.flags, &[true, false, false]);
        assert_eq!(array.equations, &[0]);
        assert_eq!(array.has_non_zero(), true);
        assert_eq!(array.value(0, 0.0), 110.0);
    }

    #[test]
    fn bc_prescribed_array_works_beam_3d() {
        #[rustfmt::skip]
        let mesh = Mesh {
            ndim: 3,
            points: vec![
                Point { id: 0, marker: 0, coords: vec![0.0, 0.0, 0.0] },
                Point { id: 1, marker: 0, coords: vec![1.0, 1.0, 1.0] },
            ],
            cells: vec![
                Cell { id: 0, attribute: 1, kind: GeoKind::Lin2, points: vec![0, 1] },
            ],
        };
        let p1 = ParamBeam::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Beam(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Dof::Ux, 1.0)
            .points(&[0], Dof::Uy, 2.0)
            .points(&[0], Dof::Uz, 3.0)
            .points(&[0], Dof::Rx, 4.0)
            .points(&[0], Dof::Ry, 5.0)
            .points(&[0], Dof::Rz, 6.0);
        let array = BcPrescribed::new(&base, &essential).unwrap();
        assert_eq!(
            array.flags,
            &[
                true, true, true, true, true, true, //        0 Ux,Uy,Uz, Rx,Ry,Rz
                false, false, false, false, false, false, //  1 Ux,Uy,Uz, Rx,Ry,Rz
            ]
        );
        assert_eq!(array.has_non_zero(), true);
        let mut eqs = array.equations.clone();
        eqs.sort();
        assert_eq!(&eqs, &[0, 1, 2, 3, 4, 5]);
    }

    #[test]
    fn bc_prescribed_array_works_mixed() {
        //                     {Ux→15}
        //    {Ux→21}          {Uy→16}
        //    {Uy→22}  {Ux→19} {Rz→17}
        //    {Pl→23}  {Uy→20} {Pl→18} {Ux→13}
        //         8------7------6._   {Uy→14}
        //         |       [3](3)|  '-.5
        //         |  [0]        |     '-._
        // {Ux→24} 9  (1)      *10  [1]    '4 {Ux→11}
        // {Uy→25} |             |  (2)  .-'  {Uy→12}
        //         |       [2](3)|   _.3'
        //         0------1------2.-'  {Ux→9}
        //     {Ux→0}  {Ux→3}  {Ux→5}  {Uy→10}
        //     {Uy→1}  {Uy→4}  {Uy→6}
        //     {Pl→2}          {Rz→7}
        //                     {Pl→8}
        //  *10 => {Ux→26, Uy→27, Rz→28}
        let mesh = Samples::qua8_tri6_lin2();
        let p1 = ParamPorousSldLiq::sample_brooks_corey_constant_elastic();
        let p2 = ParamSolid::sample_linear_elastic();
        let p3 = ParamBeam::sample();
        let base = FemBase::new(
            &mesh,
            [(1, Elem::PorousSldLiq(p1)), (2, Elem::Solid(p2)), (3, Elem::Beam(p3))],
        )
        .unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Dof::Ux, 0.0)
            .points(&[0], Dof::Uy, 1.0)
            .points(&[0], Dof::Pl, 2.0)
            .points(&[1], Dof::Ux, 3.0)
            .points(&[1], Dof::Uy, 4.0)
            .points(&[2], Dof::Ux, 5.0)
            .points(&[2], Dof::Uy, 6.0)
            .points(&[2], Dof::Rz, 7.0)
            .points(&[2], Dof::Pl, 8.0)
            .points(&[3], Dof::Ux, 9.0)
            .points(&[3], Dof::Uy, 10.0)
            .points(&[4], Dof::Ux, 11.0)
            .points(&[4], Dof::Uy, 12.0)
            .points(&[5], Dof::Ux, 13.0)
            .points(&[5], Dof::Uy, 14.0)
            .points(&[6], Dof::Ux, 15.0)
            .points(&[6], Dof::Uy, 16.0)
            .points(&[6], Dof::Rz, 17.0)
            .points(&[6], Dof::Pl, 18.0)
            .points(&[7], Dof::Ux, 19.0)
            .points(&[7], Dof::Uy, 20.0)
            .points(&[8], Dof::Ux, 21.0)
            .points(&[8], Dof::Uy, 22.0)
            .points(&[8], Dof::Pl, 23.0)
            .points(&[9], Dof::Ux, 24.0)
            .points(&[9], Dof::Uy, 25.0)
            .points(&[10], Dof::Ux, 26.0)
            .points(&[10], Dof::Uy, 27.0)
            .points(&[10], Dof::Rz, 28.0);
        let _array = BcPrescribed::new(&base, &essential).unwrap();
        #[rustfmt::skip]
        let _correct = &[            // point
             0.0,  1.0,  2.0,       //  0 (Ux, 0) (Uy, 1) (Pl,2)
             3.0,  4.0,             //  1 (Ux, 3) (Uy, 4)
             5.0,  6.0,  7.0,  8.0, //  2 (Ux, 5) (Uy, 6) (Rz,7) (Pl,8)
             9.0, 10.0,             //  3 (Ux, 9) (Uy,10)
            11.0, 12.0,             //  4 (Ux,11) (Uy,12)
            13.0, 14.0,             //  5 (Ux,13) (Uy,14)
            15.0, 16.0, 17.0, 18.0, //  6 (Ux,15) (Uy,16) (Rz,17) (Pl,18)
            19.0, 20.0,             //  7 (Ux,19) (Uy,20)
            21.0, 22.0, 23.0,       //  8 (Ux,21) (Uy,22) (Pl,23)
            24.0, 25.0,             //  9 (Ux,24) (Uy,25)
            26.0, 27.0, 28.0,       // 10 (Ux,26) (Uy,27) (Rz,28)
        ];
    }

    #[test]
    fn bc_prescribed_array_works_porous_sld_liq_gas() {
        let mesh = Samples::one_tri6();
        let p1 = ParamPorousSldLiqGas::sample_brooks_corey_constant_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::PorousSldLiqGas(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Dof::Ux, 1.0)
            .points(&[0], Dof::Uy, 2.0)
            .points(&[0], Dof::Pl, 3.0)
            .points(&[0], Dof::Pg, 4.0)
            .points(&[1], Dof::Ux, 5.0)
            .points(&[1], Dof::Uy, 6.0)
            .points(&[1], Dof::Pl, 7.0)
            .points(&[1], Dof::Pg, 8.0)
            .points(&[2], Dof::Ux, 9.0)
            .points(&[2], Dof::Uy, 10.0)
            .points(&[2], Dof::Pl, 11.0)
            .points(&[2], Dof::Pg, 12.0);
        let _array = BcPrescribed::new(&base, &essential).unwrap();
        #[rustfmt::skip]
        let _correct = &[
            1.0,  2.0,  3.0,  4.0, // 0 Ux,Uy,Pl,Pg
            5.0,  6.0,  7.0,  8.0, // 1 Ux,Uy,Pl,Pg
            9.0, 10.0, 11.0, 12.0, // 2 Ux,Uy,Pl,Pg
            0.0,  0.0,             // 3 Ux,Uy
            0.0,  0.0,             // 4 Ux,Uy
            0.0,  0.0,             // 5 Ux,Uy
        ];
    }

    #[test]
    fn bc_prescribed_array_works_triangles() {
        //       {4} 4---.__
        //          / \     `--.___3 {3}  [#] indicates id
        //         /   \          / \     (#) indicates attribute
        //        /     \  [1]   /   \    {#} indicates equation id
        //       /  [0]  \ (1)  / [2] \
        //      /   (1)   \    /  (1)  \
        // {0} 0---.__     \  /      ___2 {2}
        //            `--.__\/__.---'
        //               {1} 1
        let mesh = Samples::three_tri3();
        let p1 = ParamPorousLiq::sample_brooks_corey_constant();
        let base = FemBase::new(&mesh, [(1, Elem::PorousLiq(p1))]).unwrap();
        let mut essential = Essential::new();
        essential.points(&[0, 4], Dof::Pl, 0.0);
        let values = BcPrescribed::new(&base, &essential).unwrap();
        assert_eq!(values.flags, &[true, false, false, false, true]);
        let mut eqs = values.equations.clone();
        eqs.sort();
        assert_eq!(eqs, &[0, 4]);

        //       {8} 4---.__
        //       {9}/ \     `--.___3 {6}   [#] indicates id
        //         /   \          / \{7}   (#) indicates attribute
        //        /     \  [1]   /   \     {#} indicates equation number
        //       /  [0]  \ (1)  / [2] \
        // {0}  /   (1)   \    /  (1)  \
        // {1} 0---.__     \  /      ___2 {4}
        //            `--.__\/__.---'     {5}
        //                   1 {2}
        //                     {3}
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Dof::Ux, 0.0)
            .points(&[0], Dof::Uy, 0.0)
            .points(&[1, 2], Dof::Uy, 0.0);
        let values = BcPrescribed::new(&base, &essential).unwrap();
        assert_eq!(
            values.flags,
            //   0     1      2     3      4     5      6      7      8      9
            &[true, true, false, true, false, true, false, false, false, false]
        );
        let mut eqs = values.equations.clone();
        eqs.sort();
        assert_eq!(eqs, &[0, 1, 3, 5]);
    }
}
