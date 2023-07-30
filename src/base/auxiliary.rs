use russell_tensor::{Mandel, Tensor2};

/// Allocates a symmetric Tensor2
pub fn new_tensor2(two_dim: bool) -> Tensor2 {
    if two_dim {
        Tensor2::new(Mandel::Symmetric2D)
    } else {
        Tensor2::new(Mandel::Symmetric)
    }
}

/// Allocates a symmetric Tensor2 given the space dimension
pub fn new_tensor2_ndim(space_ndim: usize) -> Tensor2 {
    if space_ndim == 2 {
        Tensor2::new(Mandel::Symmetric2D)
    } else {
        Tensor2::new(Mandel::Symmetric)
    }
}
