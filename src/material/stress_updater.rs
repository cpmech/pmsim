#![allow(unused)]

use super::StressStrainModel;
use russell_tensor::Tensor2;

/// Simulates the stress-update process at the material point (Gauss point) level
pub struct StressUpdater {
    model: StressStrainModel,
}

impl StressUpdater {
    /// Allocates a new instance
    pub fn new(model: StressStrainModel) -> Self {
        StressUpdater { model }
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressUpdater;
    use crate::base::SampleParams;
    use crate::material::StressStrainModel;

    #[test]
    fn new_works() {
        let param = SampleParams::param_solid();
        let model = StressStrainModel::new(&param, false, false).unwrap();
        let mut updater = StressUpdater::new(model);
    }
}
