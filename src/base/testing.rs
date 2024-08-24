use super::Config;
use gemlab::mesh::Mesh;

#[allow(dead_code)]
pub(crate) fn new_empty_mesh_2d() -> Mesh {
    Mesh {
        ndim: 2,
        points: Vec::new(),
        cells: Vec::new(),
    }
}

#[allow(dead_code)]
pub(crate) fn new_empty_mesh_3d() -> Mesh {
    Mesh {
        ndim: 3,
        points: Vec::new(),
        cells: Vec::new(),
    }
}

#[allow(dead_code)]
pub(crate) fn new_empty_config_2d() -> Config {
    let mesh = new_empty_mesh_2d();
    Config::new(&mesh)
}

#[allow(dead_code)]
pub(crate) fn new_empty_config_3d() -> Config {
    let mesh = new_empty_mesh_3d();
    Config::new(&mesh)
}
