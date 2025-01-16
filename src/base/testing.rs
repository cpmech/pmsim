use gemlab::mesh::Mesh;
use russell_lab::Vector;
use russell_tensor::{Mandel, Tensor2};

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

/// Generates a displacement field corresponding to a horizontal stretching
///
/// Note: This function only works for a homogeneous mesh; with same element kinds.
///
/// # Input
///
/// * `mesh` -- the mesh
/// * `eps_xx` -- the horizontal component of strain ε_xx
#[allow(dead_code)]
pub(crate) fn generate_horizontal_displacement_field(mesh: &Mesh, eps_xx: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    for p in 0..npoint {
        let x = mesh.points[p].coords[0];
        uu[0 + mesh.ndim * p] = eps_xx * x;
    }
    uu
}

/// Generates a displacement field corresponding to a vertical stretching
///
/// Note: This function only works for a homogeneous mesh; with same element kinds.
///
/// # Input
///
/// * `mesh` -- the mesh
/// * `eps_yy` -- the vertical component of strain ε_yy
#[allow(dead_code)]
pub(crate) fn generate_vertical_displacement_field(mesh: &Mesh, eps_yy: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    for p in 0..npoint {
        let y = mesh.points[p].coords[1];
        uu[1 + mesh.ndim * p] = eps_yy * y;
    }
    uu
}

/// Generates a displacement field corresponding to a simple shear deformation
///
/// Note: This function only works for a homogeneous mesh; with same element kinds.
///
/// # Input
///
/// * `mesh` -- the mesh
/// * `eps_xy` -- the shear component of strain ε_xy
#[allow(dead_code)]
pub(crate) fn generate_shear_displacement_field(mesh: &Mesh, eps_xy: f64) -> Vector {
    let npoint = mesh.points.len();
    let mut uu = Vector::new(mesh.ndim * npoint);
    let gamma_xy = 2.0 * eps_xy; // ε_xy = 𝛾_xy / 2
    for p in 0..npoint {
        let y = mesh.points[p].coords[1];
        uu[0 + mesh.ndim * p] = gamma_xy * y;
    }
    uu
}

/// Returns the elastic solution (plane-strain or 3D) corresponding to a horizontal stretching
///
/// Returns `(strain, stress)`
///
/// # Input
///
/// * `young` -- Young's modulus
/// * `poisson` -- Poisson's coefficient
/// * `ndim` -- the space dimension
/// * `eps_xx` -- the horizontal component of strain ε_xx
#[allow(dead_code)]
pub(crate) fn elastic_solution_horizontal_displacement_field(
    young: f64,
    poisson: f64,
    ndim: usize,
    eps_xx: f64,
) -> (Tensor2, Tensor2) {
    // Mandel representation
    let mandel = if ndim == 2 {
        Mandel::Symmetric2D
    } else {
        Mandel::Symmetric
    };
    // strain tensor
    let ______ = 0.0;
    let strain = Tensor2::from_matrix(
        &[
            [eps_xx, ______, ______],
            [______, ______, ______],
            [______, ______, ______],
        ],
        mandel,
    )
    .unwrap();
    // stress tensor
    let c = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let sig_xx = c * eps_xx * (1.0 - poisson);
    let sig_yy = c * eps_xx * poisson;
    let sig_zz = c * eps_xx * poisson;
    let stress = Tensor2::from_matrix(
        &[
            [sig_xx, ______, ______],
            [______, sig_yy, ______],
            [______, ______, sig_zz],
        ],
        mandel,
    )
    .unwrap();
    (strain, stress)
}

/// Returns the elastic solution (plane-strain or 3D) corresponding to a vertical stretching
///
/// Returns `(strain, stress)`
///
/// # Input
///
/// * `young` -- Young's modulus
/// * `poisson` -- Poisson's coefficient
/// * `ndim` -- the space dimension
/// * `eps_yy` -- the vertical component of strain ε_yy
#[allow(dead_code)]
pub(crate) fn elastic_solution_vertical_displacement_field(
    young: f64,
    poisson: f64,
    ndim: usize,
    eps_yy: f64,
) -> (Tensor2, Tensor2) {
    // Mandel representation
    let mandel = if ndim == 2 {
        Mandel::Symmetric2D
    } else {
        Mandel::Symmetric
    };
    // strain tensor
    let ______ = 0.0;
    let strain = Tensor2::from_matrix(
        &[
            [______, ______, ______],
            [______, eps_yy, ______],
            [______, ______, ______],
        ],
        mandel,
    )
    .unwrap();
    // stress tensor
    let c = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let sig_xx = c * eps_yy * poisson;
    let sig_yy = c * eps_yy * (1.0 - poisson);
    let sig_zz = c * eps_yy * poisson;
    let stress = Tensor2::from_matrix(
        &[
            [sig_xx, ______, ______],
            [______, sig_yy, ______],
            [______, ______, sig_zz],
        ],
        mandel,
    )
    .unwrap();
    (strain, stress)
}

/// Returns the elastic solution (plane-strain or 3D) corresponding to a simple shear deformation
///
/// Returns `(strain, stress)`
///
/// # Input
///
/// * `young` -- Young's modulus
/// * `poisson` -- Poisson's coefficient
/// * `ndim` -- the space dimension
/// * `eps_xy` -- the shear component of strain ε_xy
#[allow(dead_code)]
pub(crate) fn elastic_solution_shear_displacement_field(
    young: f64,
    poisson: f64,
    ndim: usize,
    eps_xy: f64,
) -> (Tensor2, Tensor2) {
    // Mandel representation
    let mandel = if ndim == 2 {
        Mandel::Symmetric2D
    } else {
        Mandel::Symmetric
    };
    // strain tensor
    let ______ = 0.0;
    let strain = Tensor2::from_matrix(
        &[
            [______, eps_xy, ______],
            [eps_xy, ______, ______],
            [______, ______, ______],
        ],
        mandel,
    )
    .unwrap();
    // stress tensor
    let c = young / ((1.0 + poisson) * (1.0 - 2.0 * poisson));
    let sig_xy = c * (1.0 - 2.0 * poisson) * eps_xy;
    let stress = Tensor2::from_matrix(
        &[
            [______, sig_xy, ______],
            [sig_xy, ______, ______],
            [______, ______, ______],
        ],
        mandel,
    )
    .unwrap();
    (strain, stress)
}
