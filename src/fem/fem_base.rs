use crate::base::{Attributes, Elem, ElementDofsMap, Equations};
use crate::StrError;
use gemlab::mesh::{Cell, CellAttribute, Mesh};
use serde::{Deserialize, Serialize};
use std::ffi::OsStr;
use std::fs::{self, File};
use std::io::BufReader;
use std::path::Path;

/// Holds the material parameters, element attributes, and equation numbers
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct FemBase {
    /// Holds all attributes
    pub attributes: Attributes,

    /// Holds the element information such as local DOFs and equation numbers
    pub information: ElementDofsMap,

    /// Holds all DOF numbers
    pub equations: Equations,
}

impl FemBase {
    /// Allocates a new instance
    pub fn new<const N: usize>(mesh: &Mesh, arr: [(CellAttribute, Elem); N]) -> Result<Self, StrError> {
        let attributes = Attributes::from(arr);
        let information = ElementDofsMap::new(&mesh, &attributes)?;
        let equations = Equations::new(&mesh, &information).unwrap(); // cannot fail
        Ok(FemBase {
            attributes,
            information,
            equations,
        })
    }

    /// Returns the number of local equations
    pub fn n_local_eq(&self, cell: &Cell) -> Result<usize, StrError> {
        let info = self.information.get(cell)?;
        Ok(info.n_equation)
    }

    /// Reads a JSON file containing the base data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn read_json<P>(full_path: &P) -> Result<Self, StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        let data = File::open(path).map_err(|_| "cannot open base file")?;
        let buffered = BufReader::new(data);
        let state = serde_json::from_reader(buffered).map_err(|_| "cannot parse base file")?;
        Ok(state)
    }

    /// Writes a JSON file with the base data
    ///
    /// # Input
    ///
    /// * `full_path` -- may be a String, &str, or Path
    pub fn write_json<P>(&self, full_path: &P) -> Result<(), StrError>
    where
        P: AsRef<OsStr> + ?Sized,
    {
        let path = Path::new(full_path).to_path_buf();
        if let Some(p) = path.parent() {
            fs::create_dir_all(p).map_err(|_| "cannot create directory")?;
        }
        let mut file = File::create(&path).map_err(|_| "cannot create base file")?;
        serde_json::to_writer(&mut file, &self).map_err(|_| "cannot write base file")?;
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::FemBase;
    use crate::base::{Elem, ParamDiffusion, ParamSolid};
    use gemlab::mesh::{Cell, Samples};
    use gemlab::shapes::GeoKind;

    #[test]
    fn new_handles_errors() {
        let mesh = Samples::one_tri3();
        let p2 = ParamSolid::sample_linear_elastic();
        assert_eq!(
            FemBase::new(&mesh, [(2, Elem::Solid(p2))]).err(),
            Some("cannot find CellAttribute in Attributes map")
        );
    }

    #[test]
    fn new_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        assert_eq!(base.equations.n_equation, 6);
    }

    #[test]
    fn n_local_eq_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamDiffusion::sample();
        let base = FemBase::new(&mesh, [(1, Elem::Diffusion(p1))]).unwrap();
        assert_eq!(base.n_local_eq(&mesh.cells[0]).unwrap(), 3);

        let wrong_cell = Cell {
            id: 0,
            attribute: 1,
            kind: GeoKind::Qua4,
            points: vec![0, 1, 2, 3],
        };
        assert_eq!(
            base.n_local_eq(&wrong_cell).err(),
            Some("cannot find (CellAttribute, GeoKind) in ElementDofsMap")
        );
    }

    #[test]
    fn derive_works() {
        let mesh = Samples::one_tri3();
        let p1 = ParamSolid::sample_linear_elastic();
        let base = FemBase::new(&mesh, [(1, Elem::Solid(p1))]).unwrap();
        let clone = base.clone();
        let str_ori = format!("{:?}", clone).to_string();
        assert_eq!(format!("{:?}", clone), str_ori);
        // serialize
        let json = serde_json::to_string(&clone).unwrap();
        // deserialize
        let read: FemBase = serde_json::from_str(&json).unwrap();
        assert_eq!(format!("{:?}", read.attributes), format!("{:?}", base.attributes));
        assert_eq!(format!("{:?}", read.information), format!("{:?}", base.information));
        assert_eq!(format!("{}", read.equations), format!("{}", base.equations));
    }
}
