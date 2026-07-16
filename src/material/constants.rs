/// Tolerance to detect elastic regime
pub(crate) const F_TOL: f64 = 1e-6;

/// Indicates that the simulation should not stop
pub(crate) const KEEP_RUNNING: bool = false;

/// Number of divisions for dense output during stress-strain history recording
pub(crate) const HISTORY_N_OUT: usize = 20;

/// Tolerance to avoid negative plastic numerator (df/dσ : Dₑ : Δε)
pub(crate) const NUMERATOR_TOL: f64 = 1e-8;

/// Holds the tolerance to truncate the Chebyshev series used in root-finding
pub(crate) const CHEBYSHEV_TOL: f64 = 1e-8;

/// Holds the pseudo-time tolerance
pub(crate) const PSEUDO_TIME_TOL: f64 = 1e-7;
