use super::FemBase;
use crate::base::Essential;
use crate::StrError;
use russell_lab::Vector;

/// Assists in calculating a prescribed value boundary condition
///
/// This data structure corresponds to a single Essential (Dirichlet) boundary condition
pub struct BcPrescribed<'a> {
    /// Equation corresponding to the prescribed value
    eq: usize,

    /// Specified BC value (overridden by the function, if not None)
    value: f64,

    /// Function to calculate the BC value (overrides the value, if not None)
    function: Option<&'a Box<dyn Fn(f64) -> f64 + 'a>>,
}

/// Implements an array of BcPrescribed
pub struct BcPrescribedArray<'a> {
    /// All values
    pub all: Vec<BcPrescribed<'a>>,

    /// An array indicating which DOFs (equations) are prescribed
    ///
    /// The length of `flags` is equal to `n_equation`, the total number of DOFs (total number of equations).
    pub flags: Vec<bool>,

    /// Array with only the DOFs numbers of the prescribed equations
    ///
    /// Compared to the array `flags`, this is a "smaller" array with only the prescribed DOFs numbers.
    pub equations: Vec<usize>,
}

impl<'a> BcPrescribed<'a> {
    /// Sets prescribed value in the solution vector
    pub fn set_value(&self, duu: &mut Vector, uu: &mut Vector, time: f64) {
        let value = match self.function {
            Some(f) => (f)(time),
            None => self.value,
        };
        duu[self.eq] = value - uu[self.eq];
        uu[self.eq] = value;
    }
}

impl<'a> BcPrescribedArray<'a> {
    /// Allocates a new instance
    pub fn new(base: &FemBase, essential: &'a Essential) -> Result<Self, StrError> {
        let mut all = Vec::new();
        let mut flags = vec![false; base.equations.n_equation];
        let mut equations = Vec::new();
        for ((point_id, dof), (value, f_index)) in &essential.all {
            let eq = base.equations.eq(*point_id, *dof)?;
            let function = match f_index {
                Some(index) => Some(&essential.functions[*index]),
                None => None,
            };
            all.push(BcPrescribed {
                eq,
                value: *value,
                function,
            });
            flags[eq] = true;
            equations.push(eq);
        }
        Ok(BcPrescribedArray { all, flags, equations })
    }

    /// Sets all prescribed values in the solution vector
    pub fn apply(&self, duu: &mut Vector, uu: &mut Vector, time: f64) {
        self.all.iter().for_each(|e| e.set_value(duu, uu, time));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::BcPrescribedArray;
    use crate::base::{Dof, Elem, Essential, ParamBeam, ParamDiffusion};
    use crate::base::{ParamPorousLiq, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
    use crate::fem::FemBase;
    use gemlab::mesh::{Cell, Mesh, Point, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();

        let mut essential = Essential::new();
        essential.points(&[100], Dof::Ux, 0.0);
        assert_eq!(
            BcPrescribedArray::new(&base, &essential).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );

        let mut essential = Essential::new();
        essential.points(&[0], Dof::T, 0.0);
        assert_eq!(
            BcPrescribedArray::new(&base, &essential).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );
    }

    #[test]
    fn set_values_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        let mut essential = Essential::new();
        essential.points(&[0], Dof::T, 110.0);
        let mut duu = Vector::new(base.equations.n_equation);
        let mut uu = Vector::new(base.equations.n_equation);
        uu.fill(100.0);
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
        values.apply(&mut duu, &mut uu, 0.0);
        assert_eq!(duu.as_data(), &[10.0, 0.0, 0.0]);
        assert_eq!(uu.as_data(), &[110.0, 100.0, 100.0]);
    }

    #[test]
    fn set_values_works_beam_3d() {
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
        let mut duu = Vector::new(base.equations.n_equation);
        let mut uu = Vector::new(base.equations.n_equation);
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
        values.apply(&mut duu, &mut uu, 0.0);
        #[rustfmt::skip]
        let correct = &[
            1.0, 2.0, 3.0, 4.0, 5.0, 6.0, //  0 Ux,Uy,Uz, Rx,Ry,Rz
            0.0, 0.0, 0.0, 0.0, 0.0, 0.0, //  1 Ux,Uy,Uz, Rx,Ry,Rz
        ];
        assert_eq!(duu.as_data(), correct);
        assert_eq!(uu.as_data(), correct);
    }

    #[test]
    fn set_values_works_mixed() {
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
        let mut duu = Vector::new(base.equations.n_equation);
        let mut uu = Vector::new(base.equations.n_equation);
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
        values.apply(&mut duu, &mut uu, 0.0);
        #[rustfmt::skip]
        let correct = &[            // point
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
        assert_eq!(duu.as_data(), correct);
        assert_eq!(uu.as_data(), correct);
    }

    #[test]
    fn set_values_works_porous_sld_liq_gas() {
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
        let mut duu = Vector::new(base.equations.n_equation);
        let mut uu = Vector::new(base.equations.n_equation);
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
        values.apply(&mut duu, &mut uu, 0.0);
        #[rustfmt::skip]
        let correct = &[
            1.0,  2.0,  3.0,  4.0, // 0 Ux,Uy,Pl,Pg
            5.0,  6.0,  7.0,  8.0, // 1 Ux,Uy,Pl,Pg
            9.0, 10.0, 11.0, 12.0, // 2 Ux,Uy,Pl,Pg
            0.0,  0.0,             // 3 Ux,Uy
            0.0,  0.0,             // 4 Ux,Uy
            0.0,  0.0,             // 5 Ux,Uy
        ];
        assert_eq!(duu.as_data(), correct);
        assert_eq!(uu.as_data(), correct);
    }

    #[test]
    fn prescribed_arrays_are_correct() {
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
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
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
        let values = BcPrescribedArray::new(&base, &essential).unwrap();
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
