use std::path::{Path, PathBuf};

pub struct FilePath {}

impl FilePath {
    /// Returns the filepath of a mesh file (.mesh binary format)
    ///
    /// # Input
    ///
    /// * `filename_key` -- the filename without path and extension; ".mesh" will added
    /// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/meshes" directory
    pub fn mesh(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
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
    pub fn msh(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
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
    pub fn svg(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
        let mut filename = String::from(filename_key);
        filename.push_str(".svg");
        if use_tmp_dir {
            Path::new("/tmp/pmsim").join(filename)
        } else {
            Path::new("data").join("figures").join(filename)
        }
    }

    /// Returns the filepath of a figure (.svg) file
    ///
    /// # Input
    ///
    /// * `filename_key` -- the filename without path and extension; ".svg" will added
    /// * `suffix` -- adds an extra word to the filename, e.g., `_mesh`
    /// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/figures" directory
    pub fn svg_suffix(filename_key: &str, suffix: &str, use_tmp_dir: bool) -> PathBuf {
        let mut filename = String::from(filename_key);
        filename.push_str(suffix);
        filename.push_str(".svg");
        if use_tmp_dir {
            Path::new("/tmp/pmsim").join(filename)
        } else {
            Path::new("data").join("figures").join(filename)
        }
    }
}
