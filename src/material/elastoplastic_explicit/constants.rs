/// Indicates that the simulation should not stop
pub(super) const KEEP_RUNNING: bool = false;

/// Number of divisions for dense output during stress-strain history recording
pub(super) const HISTORY_N_OUT: usize = 20;

/// Tolerance to avoid negative plastic numerator (df/dσ : Dₑ : Δε)
pub(super) const NUMERATOR_TOL: f64 = 1e-8;

/// Holds the tolerance to truncate the Chebyshev series used in root-finding
pub(super) const CHEBYSHEV_TOL: f64 = 1e-8;

/// Holds the pseudo-time tolerance
pub(super) const PSEUDO_TIME_TOL: f64 = 1e-7;
