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

    /// Returns the filepath of a Paraview (.vtu) file
    ///
    /// # Input
    ///
    /// * `filename_key` -- the filename without path and extension; ".vtu" will added
    /// * `use_tmp_dir` -- use "/tmp/pmsim" instead of local "data/paraview" directory
    pub fn vtu(filename_key: &str, use_tmp_dir: bool) -> PathBuf {
        let mut filename = String::from(filename_key);
        filename.push_str(".vtu");
        if use_tmp_dir {
            Path::new("/tmp/pmsim").join(filename)
        } else {
            Path::new("data").join("paraview").join(filename)
        }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FilePath;
    use std::ffi::OsStr;

    #[test]
    fn paths_are_correct() {
        assert_eq!(
            FilePath::mesh("the_cube", false).as_os_str(),
            OsStr::new("data/meshes/the_cube.mesh")
        );
        assert_eq!(
            FilePath::mesh("the_cube", true).as_os_str(),
            OsStr::new("/tmp/pmsim/the_cube.mesh")
        );

        assert_eq!(
            FilePath::msh("the_cube", false).as_os_str(),
            OsStr::new("data/meshes/the_cube.msh")
        );
        assert_eq!(
            FilePath::msh("the_cube", true).as_os_str(),
            OsStr::new("/tmp/pmsim/the_cube.msh")
        );

        assert_eq!(
            FilePath::svg("the_cube", false).as_os_str(),
            OsStr::new("data/figures/the_cube.svg")
        );
        assert_eq!(
            FilePath::svg("the_cube", true).as_os_str(),
            OsStr::new("/tmp/pmsim/the_cube.svg")
        );

        assert_eq!(
            FilePath::svg_suffix("the_cube", "_mesh", false).as_os_str(),
            OsStr::new("data/figures/the_cube_mesh.svg")
        );
        assert_eq!(
            FilePath::svg_suffix("the_cube", "_mesh", true).as_os_str(),
            OsStr::new("/tmp/pmsim/the_cube_mesh.svg")
        );

        assert_eq!(
            FilePath::vtu("the_cube", false).as_os_str(),
            OsStr::new("data/paraview/the_cube.vtu")
        );
        assert_eq!(
            FilePath::vtu("the_cube", true).as_os_str(),
            OsStr::new("/tmp/pmsim/the_cube.vtu")
        );
    }
}
