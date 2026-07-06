//! This module contains analytical solutions to some problems or reference solutions for testing and verifications

mod confined_flow;
mod elast_plane_strain_flexible_foot;
mod elast_plane_strain_pres_cylin;
mod plast_circular_plate;
mod polar_coordinates;
mod pres_cylin_plane_strain;
mod pres_sphere_axisymmetric;

pub use confined_flow::*;
pub use elast_plane_strain_flexible_foot::*;
pub use elast_plane_strain_pres_cylin::*;
pub use plast_circular_plate::*;
pub use polar_coordinates::*;
pub use pres_cylin_plane_strain::*;
pub use pres_sphere_axisymmetric::*;
