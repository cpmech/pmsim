/// Defines the directory where the simulation result files are saved
pub const DEFAULT_OUT_DIR: &str = "/tmp/pmsim/results";

/// Defines an auxiliary directory where the test result files are saved
pub const DEFAULT_TEST_DIR: &str = "/tmp/pmsim/test";

/// Holds the number of internal values for the LinearElastic model
///
/// `internal_values = []`
pub const N_INT_VAL_LINEAR_ELASTIC: usize = 0;

/// Holds the number of internal values for the VonMises model
///
/// `internal_values = [z, lambda]`
///
/// * `z` -- is the size of yield surface
/// * `lambda` -- is the algorithmic Lagrange multiplier
pub const N_INT_VAL_VON_MISES: usize = 2;

/// Holds the number of internal values for the DruckerPrager model
///
/// `internal_values = [z, lambda, apex_return]`
///
/// * `z` -- is the size of yield surface
/// * `lambda` -- is the algorithmic Lagrange multiplier
/// * `apex_return` -- indicates that the return algorithm handled the apex discontinuity directly
pub const N_INT_VAL_DRUCKER_PRAGER: usize = 3;

/// Holds the number of internal values for the CamClay model
///
/// `internal_values = [z]`
///
/// * `z` -- is the size of yield surface
pub const N_INT_VAL_CAM_CLAY: usize = 1;
