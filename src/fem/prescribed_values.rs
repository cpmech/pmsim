use super::Data;
use crate::base::{Ebc, Essential};
use crate::StrError;
use gemlab::mesh::{Point, PointId};
use russell_lab::Vector;

/// Assists in calculating prescribed values
pub struct PrescribedValue<'a> {
    pub point: &'a Point,
    pub ebc: Ebc,
    pub eq: usize,
}

/// Holds a collection of prescribed (primary) values
pub struct PrescribedValues<'a> {
    pub all: Vec<PrescribedValue<'a>>,
}

impl<'a> PrescribedValue<'a> {
    /// Allocates new instance
    pub fn new(data: &'a Data, point_id: PointId, ebc: Ebc) -> Result<Self, StrError> {
        if point_id >= data.mesh.points.len() {
            return Err("cannot initialize prescribed value because PointId is out-of-bounds");
        }
        Ok(PrescribedValue {
            point: &data.mesh.points[point_id],
            ebc,
            eq: data.equations.eq(point_id, ebc.dof())?,
        })
    }

    /// Sets prescribed value in the solution vector
    pub fn set_value(&self, uu: &mut Vector, time: f64) {
        let value = match self.ebc {
            Ebc::Ux(f) => f(time),
            Ebc::Uy(f) => f(time),
            Ebc::Uz(f) => f(time),
            Ebc::Rx(f) => f(time),
            Ebc::Ry(f) => f(time),
            Ebc::Rz(f) => f(time),
            Ebc::T(f) => f(time),
            Ebc::Pl(f) => f(time),
            Ebc::Pg(f) => f(time),
        };
        uu[self.eq] = value;
    }
}

impl<'a> PrescribedValues<'a> {
    /// Allocates new instance
    pub fn new(data: &'a Data, bcs: &Essential) -> Result<Self, StrError> {
        let res: Result<Vec<_>, _> = bcs
            .all
            .iter()
            .map(|((point_id, _), ebc)| PrescribedValue::new(data, *point_id, *ebc))
            .collect();
        match res {
            Ok(all) => Ok(PrescribedValues { all }),
            Err(e) => Err(e),
        }
    }

    /// Sets all prescribed values in the solution vector
    #[inline]
    pub fn set_values(&self, uu: &mut Vector, time: f64) {
        self.all.iter().for_each(|e| e.set_value(uu, time));
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{PrescribedValue, PrescribedValues};
    use crate::base::{Ebc, Element, Essential, SampleParams};
    use crate::fem::Data;
    use gemlab::mesh::Samples;
    use russell_lab::Vector;

    #[test]
    fn new_captures_errors() {
        let mesh = Samples::one_tri3();
        let p1 = SampleParams::param_solid();
        let data = Data::new(&mesh, [(1, Element::Solid(p1))]).unwrap();
        let f = |_| 123.0;
        assert_eq!(f(0.0), 123.0);
        assert_eq!(
            PrescribedValue::new(&data, 123, Ebc::Ux(f)).err(),
            Some("cannot initialize prescribed value because PointId is out-of-bounds")
        );
        assert_eq!(
            PrescribedValue::new(&data, 0, Ebc::T(f)).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );

        let mut essential = Essential::new();
        essential.at(&[100], Ebc::Ux(f));
        assert_eq!(
            PrescribedValues::new(&data, &essential).err(),
            Some("cannot initialize prescribed value because PointId is out-of-bounds")
        );
        let mut essential = Essential::new();
        essential.at(&[0], Ebc::T(f));
        assert_eq!(
            PrescribedValues::new(&data, &essential).err(),
            Some("cannot find equation number corresponding to (PointId,DOF)")
        );
    }

    #[test]
    fn set_values_works() {
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
        let p1 = SampleParams::param_porous_sld_liq();
        let p2 = SampleParams::param_solid();
        let p3 = SampleParams::param_beam();
        let data = Data::new(
            &mesh,
            [
                (1, Element::PorousSldLiq(p1)),
                (2, Element::Solid(p2)),
                (3, Element::Beam(p3)),
            ],
        )
        .unwrap();
        let mut essential = Essential::new();
        essential
            .at(&[0], Ebc::Ux(|_| 0.0))
            .at(&[0], Ebc::Uy(|_| 1.0))
            .at(&[0], Ebc::Pl(|_| 2.0))
            .at(&[1], Ebc::Ux(|_| 3.0))
            .at(&[1], Ebc::Uy(|_| 4.0))
            .at(&[2], Ebc::Ux(|_| 5.0))
            .at(&[2], Ebc::Uy(|_| 6.0))
            .at(&[2], Ebc::Rz(|_| 7.0))
            .at(&[2], Ebc::Pl(|_| 8.0))
            .at(&[3], Ebc::Ux(|_| 9.0))
            .at(&[3], Ebc::Uy(|_| 10.0))
            .at(&[4], Ebc::Ux(|_| 11.0))
            .at(&[4], Ebc::Uy(|_| 12.0))
            .at(&[5], Ebc::Ux(|_| 13.0))
            .at(&[5], Ebc::Uy(|_| 14.0))
            .at(&[6], Ebc::Ux(|_| 15.0))
            .at(&[6], Ebc::Uy(|_| 16.0))
            .at(&[6], Ebc::Rz(|_| 17.0))
            .at(&[6], Ebc::Pl(|_| 18.0))
            .at(&[7], Ebc::Ux(|_| 19.0))
            .at(&[7], Ebc::Uy(|_| 20.0))
            .at(&[8], Ebc::Ux(|_| 21.0))
            .at(&[8], Ebc::Uy(|_| 22.0))
            .at(&[8], Ebc::Pl(|_| 23.0))
            .at(&[9], Ebc::Ux(|_| 24.0))
            .at(&[9], Ebc::Uy(|_| 25.0))
            .at(&[10], Ebc::Ux(|_| 26.0))
            .at(&[10], Ebc::Uy(|_| 27.0))
            .at(&[10], Ebc::Rz(|_| 28.0));
        let mut uu = Vector::new(data.equations.n_equation);
        let values = PrescribedValues::new(&data, &essential).unwrap();
        values.set_values(&mut uu, 0.0);
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
        assert_eq!(uu.as_data(), correct);
    }
}
