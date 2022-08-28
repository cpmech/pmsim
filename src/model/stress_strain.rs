use crate::StrError;
use russell_tensor::{Tensor2, Tensor4};

pub trait StressStrain {
    fn stiffness(&self, dd: &mut Tensor4, sigma: &Tensor2) -> Result<(), StrError>;
}
