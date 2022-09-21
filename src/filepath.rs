use std::path::{Path, PathBuf};

/// Returns the filepath of a mesh file (.mesh binary format)
///
/// # Input
///
/// * `filename_key` -- the filename without path and extension; ".mesh" will added
/// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/meshes" directory
pub fn filepath_mesh(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
    let mut filename = String::from(filename_key);
    filename.push_str(".mesh");
    if use_tmp_dir {
        Path::new("/tmp/pmsim").join(filename)
    } else {
        Path::new("data").join("meshes").join(filename)
    }
}

/// Returns the filepath of a mesh file (.msh text format)
///
/// # Input
///
/// * `filename_key` -- the filename without path and extension; ".msh" will added
/// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/meshes" directory
pub fn filepath_msh(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
    let mut filename = String::from(filename_key);
    filename.push_str(".msh");
    if use_tmp_dir {
        Path::new("/tmp/pmsim").join(filename)
    } else {
        Path::new("data").join("meshes").join(filename)
    }
}

/// Returns the filepath of a figure (.svg) file
///
/// # Input
///
/// * `filename_key` -- the filename without path and extension; ".svg" will added
/// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/figures" directory
pub fn filepath_svg(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
    let mut filename = String::from(filename_key);
    filename.push_str(".svg");
    if use_tmp_dir {
        Path::new("/tmp/pmsim").join(filename)
    } else {
        Path::new("data").join("figures").join(filename)
    }
}
