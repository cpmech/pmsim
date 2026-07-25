/// Holds the number of internal variables for the LinearElastic model
///
/// None
pub const NZ_LINEAR_ELASTIC: usize = 0;

/// Holds the number of internal variables for the VonMises model
///
/// The set of internal variables is z = {κ, α} where:
/// * κ is the size of the yield surface
/// * α is the accumulated norm of the deviatoric plastic strain
pub const NZ_VON_MISES: usize = 2;

/// Holds the number of internal variables for the VonMises (with Softening) model
///
/// The set of internal variables is z = {κ, α} where:
/// * κ is the size of the yield surface
/// * α is the accumulated norm of the deviatoric plastic strain
pub const NZ_VON_MISES_SOFT: usize = 2;

/// Holds the number of internal variables for the DruckerPrager model
///
/// The set of internal variables is z = {κ, α} where:
/// * κ is the size of the yield surface
/// * α is the accumulated mean plastic strain
pub const NZ_DRUCKER_PRAGER: usize = 2;

/// Holds the number of internal variables for the CamClay model
///
/// The set of internal variables is z = {κ} where:
/// * κ is the size of the yield surface
pub const NZ_CAM_CLAY: usize = 1;
