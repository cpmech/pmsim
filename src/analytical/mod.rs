//! This module contains analytical solutions to some problems or reference solutions for testing and verifications

mod elast_plane_strain_flexible_foot;
mod elast_plane_strain_pres_cylin;
mod plast_plane_strain_pres_cylin;
mod polar_coordinates;

pub use elast_plane_strain_flexible_foot::*;
pub use elast_plane_strain_pres_cylin::*;
pub use plast_plane_strain_pres_cylin::*;
pub use polar_coordinates::*;
