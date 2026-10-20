/// Holds the number of internal variables for the LinearElastic model
///
/// None
pub const NZ_LINEAR_ELASTIC: usize = 0;

/// Holds the number of internal variables for the VonMises model
///
/// The set of internal variables is z = {a, b} where:
/// * a is the accumulated plastic strain, given by ȧ = λ (isomorphic convention)
/// * b is the size of the yield surface
pub const NZ_VON_MISES: usize = 2;

/// Holds the number of internal variables for the VonMises (with Softening) model
///
/// The set of internal variables is z = {a, b} where:
/// * a is the accumulated plastic strain, given by ȧ = λ (isomorphic convention)
/// * b is the size of the yield surface
pub const NZ_VON_MISES_SOFT: usize = 2;

/// Holds the number of internal variables of the VonMisesKinematic model *besides* the back stress
///
/// The set of internal variables is z = {a, b, c} where:
/// * a is the accumulated plastic strain, given by ȧ = λ (isomorphic convention)
/// * b is the size of the yield surface
/// * c is the back stress (one tensor, thus `ncp` components)
///
/// Therefore, this model has `NZ_VON_MISES_KINEMATIC_EXTRA + ncp` internal variables.
pub const NZ_VON_MISES_KINEMATIC_EXTRA: usize = 2;

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
