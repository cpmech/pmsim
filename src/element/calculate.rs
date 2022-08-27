// use gemlab::mesh::{CellId, Mesh};
use crate::StrError;

pub trait Calculate {
    // fn new(mesh: &Mesh, cell_id: CellId) -> Self;
    fn residual(&mut self) -> Result<(), StrError>;
    fn jacobian(&mut self) -> Result<(), StrError>;
}
