/// Indicates the yield surface crossing case
#[derive(Clone, Copy, Debug)]
pub enum Case {
    /// Elastic
    AE,

    /// Elastic-elastoplastic; holds t_intersection
    AXB(f64),

    /// Elastic; going inside (with eventual crossing)
    BE,

    /// Elastic-elastoplastic; going inside then outside with two crossings; holds t_intersection
    BXP(f64),

    /// Elastoplastic
    BP,
}
