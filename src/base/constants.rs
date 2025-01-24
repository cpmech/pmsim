/// Defines the directory where the simulation result files are saved
pub const DEFAULT_OUT_DIR: &str = "/tmp/pmsim/results";

/// Defines an auxiliary directory where the test result files are saved
pub const DEFAULT_TEST_DIR: &str = "/tmp/pmsim/test";

/// Holds the number of internal variables for the LinearElastic model
///
/// `int_vars = []`
pub const NZ_LINEAR_ELASTIC: usize = 0;

/// Holds the number of internal variables for the VonMises model
///
/// `int_vars = [z, lambda]`
///
/// * `z` -- is the size of yield surface
/// * `lambda` -- is the algorithmic Lagrange multiplier
pub const NZ_VON_MISES: usize = 2;

/// Holds the number of internal variables for the DruckerPrager model
///
/// `int_vars = [z, lambda, apex_return]`
///
/// * `z` -- is the size of yield surface
/// * `lambda` -- is the algorithmic Lagrange multiplier
/// * `apex_return` -- indicates that the return algorithm handled the apex discontinuity directly
pub const NZ_DRUCKER_PRAGER: usize = 3;

/// Holds the number of internal variables for the CamClay model
///
/// `int_vars = [z]`
///
/// * `z` -- is the size of yield surface
pub const NZ_CAM_CLAY: usize = 1;
