#![allow(unused)]

use super::StressStrainModel;
use russell_tensor::Tensor2;

pub struct StressUpdater {
    model: Box<dyn StressStrainModel>,
    stress_path: Vec<Option<Tensor2>>,
    strain_path: Vec<Option<Tensor2>>,
}

impl StressUpdater {
    pub fn new(model: Box<dyn StressStrainModel>) -> Self {
        StressUpdater {
            model,
            stress_path: Vec::new(),
            strain_path: Vec::new(),
        }
    }

    pub fn push_isotropic(&mut self, sigma_m: f64, two_dim: bool) -> &mut Self {
        let mut tt = Tensor2::new_sym(two_dim);
        tt.vec[0] = sigma_m;
        tt.vec[1] = sigma_m;
        tt.vec[2] = sigma_m;
        self.stress_path.push(Some(tt));
        self.strain_path.push(None);
        self
    }

    pub fn push_octahedral(&mut self, sigma_m: f64, sigma_d: f64, lode: f64, two_dim: bool) -> &mut Self {
        let mut tt = Tensor2::new_from_oct_invariants(sigma_m, sigma_d, lode, two_dim).unwrap();
        self.stress_path.push(Some(tt));
        self.strain_path.push(None);
        self
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::StressUpdater;
    use crate::base::SampleParams;
    use crate::material::allocate_stress_strain_model;

    #[test]
    fn new_works() {
        let param = SampleParams::param_solid();
        let model = allocate_stress_strain_model(&param, false, false).unwrap();
        let mut updater = StressUpdater::new(model);

        const TWO_DIM: bool = true;

        updater
            .push_isotropic(10.0, TWO_DIM)
            .push_isotropic(20.0, TWO_DIM)
            .push_octahedral(30.0, 90.0, 1.0, TWO_DIM);
    }
}
