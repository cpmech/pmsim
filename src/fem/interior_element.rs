use super::LocalEquations;

pub struct InteriorElement {
    pub actual: Box<dyn LocalEquations>,
}
