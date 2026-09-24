/// Specifies the 2D space dimension
///
/// Note: The number recorded here is actually `DIM*2` because it corresponds to the number of
/// tensor components in the Kelvin-Mandel (**KM4**) notation for generalized planar tensors.
pub const D2: usize = 4;

/// Specifies the 3D space dimension
///
/// Note: The number recorded here is actually `DIM*2` because it corresponds to the number of
/// tensor components in the Kelvin-Mandel (**KM6**) notation for (minor-) symmetric tensors.
pub const D3: usize = 6;

/// Defines a type alias for the error type as a static string
pub type StrError = &'static str;

pub mod analytical;
pub mod base;
pub mod fem;
pub mod material;
pub mod prelude;
pub mod util;
