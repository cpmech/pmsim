use super::Config;
use gemlab::mesh::Mesh;
use russell_lab::Vector;

/// Returns a new empty 2D mesh
#[allow(dead_code)]
pub(crate) fn new_empty_mesh_2d() -> Mesh {
    Mesh {
        ndim: 2,
        points: Vec::new(),
        cells: Vec::new(),
    }
}

/// Returns a new empty 3D mesh
#[allow(dead_code)]
pub(crate) fn new_empty_mesh_3d() -> Mesh {
    Mesh {
        ndim: 3,
        points: Vec::new(),
        cells: Vec::new(),
    }
}

/// Returns a new empty configuration for 2D simulations
#[allow(dead_code)]
pub(crate) fn new_empty_config_2d() -> Config {
    let mesh = new_empty_mesh_2d();
    Config::new(&mesh)
}

/// Returns a new empty configuration for 3D simulations
#[allow(dead_code)]
pub(crate) fn new_empty_config_3d() -> Config {
    let mesh = new_empty_mesh_3d();
    Config::new(&mesh)
}

/// Generates a displacement field corresponding to a (confined) horizontal stretching
/// (only works for a homogeneous mesh; with same element kinds)
#[allow(dead_code)]
pub(crate) fn generate_horizontal_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    for p in 0..npoint {
        let x = mesh.points[p].coords[0];
        uu[0 + mesh.ndim * p] = strain * x;
    }
    uu
}

/// Generates a displacement field corresponding to a (confined) vertical stretching
/// (only works for a homogeneous mesh; with same element kinds)
#[allow(dead_code)]
pub(crate) fn generate_vertical_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    for p in 0..npoint {
        let y = mesh.points[p].coords[1];
        uu[1 + mesh.ndim * p] = strain * y;
    }
    uu
}

/// Generates a displacement field corresponding to a simple shear deformation
/// (only works for a homogeneous mesh; with same element kinds)
/// Here, strain is 𝛾; thus ε = 𝛾/2 = strain/2
#[allow(dead_code)]
pub(crate) fn generate_shear_displacement_field(mesh: &Mesh, strain: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    for p in 0..npoint {
        let y = mesh.points[p].coords[1];
        uu[0 + mesh.ndim * p] = strain * y;
    }
    uu
}
