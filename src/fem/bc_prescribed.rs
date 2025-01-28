use super::FemInput;
use crate::base::{Ebc, Essential};
use crate::StrError;
use gemlab::mesh::{Point, PointId};
use russell_lab::Vector;

/// Assists in calculating a prescribed value boundary condition
///
/// This data structure corresponds to a single Essential (Dirichlet) boundary condition
pub struct BcPrescribed<'a> {
    /// Point corresponding to the prescribed value
    pub point: &'a Point,

    /// Essential boundary condition
    pub ebc: Ebc,

    /// Equation corresponding to the prescribed value
    pub eq: usize,
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
    /// Allocates new instance
    pub fn new(input: &'a FemInput, point_id: PointId, ebc: Ebc) -> Result<Self, StrError> {
        if point_id >= input.mesh.points.len() {
            return Err("cannot initialize prescribed value because PointId is out-of-bounds");
        }
        Ok(BcPrescribed {
            point: &input.mesh.points[point_id],
            ebc,
            eq: input.equations.eq(point_id, ebc.dof())?,
        })
    }

    /// Sets prescribed value in the solution vector
    pub fn set_value(&self, duu: &mut Vector, uu: &mut Vector, _time: f64) {
        let value = match self.ebc {
            Ebc::Ux(v) => v,
            Ebc::Uy(v) => v,
            Ebc::Uz(v) => v,
            Ebc::Rx(v) => v,
            Ebc::Ry(v) => v,
            Ebc::Rz(v) => v,
            Ebc::T(v) => v,
            Ebc::Pl(v) => v,
            Ebc::Pg(v) => v,
        };
        duu[self.eq] = value - uu[self.eq];
        uu[self.eq] = value;
    }
}

impl<'a> BcPrescribedArray<'a> {
    /// Allocates new instance
    pub fn new(input: &'a FemInput, essential: &Essential) -> Result<Self, StrError> {
        let mut all = Vec::new();
        let mut flags = vec![false; input.equations.n_equation];
        let mut equations = Vec::new();
        for ((point_id, dof), ebc) in &essential.all {
            let eq = input.equations.eq(*point_id, *dof)?;
            all.push(BcPrescribed::new(input, *point_id, *ebc).unwrap()); // already checked
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
    use super::{BcPrescribed, BcPrescribedArray};
    use crate::base::{Ebc, Essential, Etype, ParamBeam, ParamDiffusion};
    use crate::base::{ParamPorousLiq, ParamPorousSldLiq, ParamPorousSldLiqGas, ParamSolid};
    use crate::fem::FemInput;
    use gemlab::mesh::{Cell, Mesh, Point, Samples};
    use gemlab::shapes::GeoKind;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))]).unwrap();
        assert_eq!(
            BcPrescribed::new(&input, 123, Ebc::Ux(0.0)).err(),
            Some("cannot initialize prescribed value because PointId is out-of-bounds")
        );
        assert_eq!(
            BcPrescribed::new(&input, 0, Ebc::T(0.0)).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );

        let mut essential = Essential::new();
        essential.points(&[100], Ebc::Ux(0.0));
        assert_eq!(
            BcPrescribedArray::new(&input, &essential).err(),
            Some("cannot find equation number because PointId is out-of-bounds")
        );
        let mut essential = Essential::new();
        essential.points(&[0], Ebc::T(0.0));
        assert_eq!(
            BcPrescribedArray::new(&input, &essential).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );
    }

    #[test]
    fn set_values_works_diffusion() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let input = FemInput::new(&mesh, [(1, Etype::Diffusion(p1))]).unwrap();
        let mut essential = Essential::new();
        essential.points(&[0], Ebc::T(110.0));
        let mut duu = Vector::new(input.equations.n_equation);
        let mut uu = Vector::new(input.equations.n_equation);
        uu.fill(100.0);
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
        let input = FemInput::new(&mesh, [(1, Etype::Beam(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Ebc::Ux(1.0))
            .points(&[0], Ebc::Uy(2.0))
            .points(&[0], Ebc::Uz(3.0))
            .points(&[0], Ebc::Rx(4.0))
            .points(&[0], Ebc::Ry(5.0))
            .points(&[0], Ebc::Rz(6.0));
        let mut duu = Vector::new(input.equations.n_equation);
        let mut uu = Vector::new(input.equations.n_equation);
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
        let input = FemInput::new(
            &mesh,
            [
                (1, Etype::PorousSldLiq(p1)),
                (2, Etype::Solid(p2)),
                (3, Etype::Beam(p3)),
            ],
        )
        .unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Ebc::Ux(0.0))
            .points(&[0], Ebc::Uy(1.0))
            .points(&[0], Ebc::Pl(2.0))
            .points(&[1], Ebc::Ux(3.0))
            .points(&[1], Ebc::Uy(4.0))
            .points(&[2], Ebc::Ux(5.0))
            .points(&[2], Ebc::Uy(6.0))
            .points(&[2], Ebc::Rz(7.0))
            .points(&[2], Ebc::Pl(8.0))
            .points(&[3], Ebc::Ux(9.0))
            .points(&[3], Ebc::Uy(10.0))
            .points(&[4], Ebc::Ux(11.0))
            .points(&[4], Ebc::Uy(12.0))
            .points(&[5], Ebc::Ux(13.0))
            .points(&[5], Ebc::Uy(14.0))
            .points(&[6], Ebc::Ux(15.0))
            .points(&[6], Ebc::Uy(16.0))
            .points(&[6], Ebc::Rz(17.0))
            .points(&[6], Ebc::Pl(18.0))
            .points(&[7], Ebc::Ux(19.0))
            .points(&[7], Ebc::Uy(20.0))
            .points(&[8], Ebc::Ux(21.0))
            .points(&[8], Ebc::Uy(22.0))
            .points(&[8], Ebc::Pl(23.0))
            .points(&[9], Ebc::Ux(24.0))
            .points(&[9], Ebc::Uy(25.0))
            .points(&[10], Ebc::Ux(26.0))
            .points(&[10], Ebc::Uy(27.0))
            .points(&[10], Ebc::Rz(28.0));
        let mut duu = Vector::new(input.equations.n_equation);
        let mut uu = Vector::new(input.equations.n_equation);
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
        let input = FemInput::new(&mesh, [(1, Etype::PorousSldLiqGas(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Ebc::Ux(1.0))
            .points(&[0], Ebc::Uy(2.0))
            .points(&[0], Ebc::Pl(3.0))
            .points(&[0], Ebc::Pg(4.0))
            .points(&[1], Ebc::Ux(5.0))
            .points(&[1], Ebc::Uy(6.0))
            .points(&[1], Ebc::Pl(7.0))
            .points(&[1], Ebc::Pg(8.0))
            .points(&[2], Ebc::Ux(9.0))
            .points(&[2], Ebc::Uy(10.0))
            .points(&[2], Ebc::Pl(11.0))
            .points(&[2], Ebc::Pg(12.0));
        let mut duu = Vector::new(input.equations.n_equation);
        let mut uu = Vector::new(input.equations.n_equation);
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
        let input = FemInput::new(&mesh, [(1, Etype::PorousLiq(p1))]).unwrap();
        let mut essential = Essential::new();
        essential.points(&[0, 4], Ebc::Pl(0.0));
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
        let input = FemInput::new(&mesh, [(1, Etype::Solid(p1))]).unwrap();
        let mut essential = Essential::new();
        essential
            .points(&[0], Ebc::Ux(0.0))
            .points(&[0], Ebc::Uy(0.0))
            .points(&[1, 2], Ebc::Uy(0.0));
        let values = BcPrescribedArray::new(&input, &essential).unwrap();
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
