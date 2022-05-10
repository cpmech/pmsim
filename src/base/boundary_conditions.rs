/// Defines natural boundary conditions (NBC)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum Nbc {
    /// Normal distributed load
    Qn,

    /// Distributed load parallel to x
    Qx,

    /// Distributed load parallel to y
    Qy,

    /// Distributed load parallel to z
    Qz,
}

/// Defines boundary conditions at points (e.g., point loads)
#[derive(Clone, Copy, Debug, Eq, PartialEq, Hash)]
pub enum BcPoint {
    /// Concentrated load parallel to x
    Fx,

    /// Concentrated load parallel to y
    Fy,

    /// Concentrated load parallel to z
    Fz,
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::{BcPoint, Nbc};

    #[test]
    fn derive_works() {
        let qn = Nbc::Qn;
        let qn_clone = qn.clone();
        assert_eq!(format!("{:?}", qn), "Qn");
        assert_eq!(qn, qn_clone);

        let fx = BcPoint::Fx;
        let fx_clone = fx.clone();
        assert_eq!(format!("{:?}", fx), "Fx");
        assert_eq!(fx, fx_clone);
    }
}
