// --- Linear Elastic Model ---

/// Holds the number of main (z) internal variables for the LinearElastic model
///
/// None
pub const NZ_LINEAR_ELASTIC: usize = 0;

/// Holds the number of extra (x) internal variables for the LinearElastic model
///
/// None
pub const NX_LINEAR_ELASTIC: usize = 0;

// --- Von Mises Model ---

/// Holds the number of main (z) internal variables for the VonMises model
///
/// `z` is the "size" of tye yield surfrace
pub const NZ_VON_MISES: usize = 1;

/// Holds the number of extra (x) internal variables for the VonMises model
///
/// `eps_bar_p` is the accumulated plastic strain (for post-processing only)
pub const NX_VON_MISES: usize = 1;

// --- Von Mises with Softening Model ---

/// Holds the number of main (z) internal variables for the VonMises (with Softening) model
///
/// `z` is the "size" of tye yield surfrace
pub const NZ_VON_MISES_SOFT: usize = 1;

/// Holds the number of extra (x) internal variables for the VonMises (with Softening) model
///
/// `eps_bar_p` is the accumulated plastic strain (needed for the softening law)
pub const NX_VON_MISES_SOFT: usize = 1;

// --- Drucker-Prager Model ---

/// Holds the number of main (z) internal variables for the DruckerPrager model
///
/// `z` is the "size" of tye yield surfrace
pub const NZ_DRUCKER_PRAGER: usize = 1;

/// Holds the number of extra (x) internal variables for the DruckerPrager model
///
/// None
pub const NX_DRUCKER_PRAGER: usize = 0;

// --- Cam-Clay Model ---

/// Holds the number of main (z) internal variables for the CamClay model
///
/// `z` is the "size" of tye yield surfrace
pub const NZ_CAM_CLAY: usize = 1;

/// Holds the number of extra (x) internal variables for the CamClay model
///
/// None
pub const NX_CAM_CLAY: usize = 0;
