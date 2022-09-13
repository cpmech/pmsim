use crate::base::ParamConductivity;
use crate::StrError;
use russell_tensor::Tensor2;

pub trait ConductivityModel: Send + Sync {
    /// Computes the conductivity tensor
    ///
    /// * `phi` may be the temperature or liquid/gas pressure
    fn tensor(&mut self, kk: &mut Tensor2, phi: f64) -> Result<(), StrError>;
}

/// Allocates conductivity model
pub fn allocate_conductivity_model(param: &ParamConductivity, two_dim: bool) -> Box<dyn ConductivityModel> {
    let model: Box<dyn ConductivityModel> = match param {
        ParamConductivity::Constant { kx, ky, kz } => Box::new(ConstantConductivityModel::new(*kx, *ky, *kz, two_dim)),
        ParamConductivity::IsotropicLinear { kr, beta } => {
            Box::new(IsotropicLinearConductivityModel::new(*kr, *beta, two_dim))
        }
        _ => panic!("todo"),
    };
    model
}

// ConstantConductivityModel  //////////////////////////////////////////////////////////////////////////////////////////

pub struct ConstantConductivityModel {
    kx: f64,
    ky: f64,
    kz: f64,
    two_dim: bool,
}

impl ConstantConductivityModel {
    pub fn new(kx: f64, ky: f64, kz: f64, two_dim: bool) -> Self {
        ConstantConductivityModel { kx, ky, kz, two_dim }
    }
}

impl ConductivityModel for ConstantConductivityModel {
    fn tensor(&mut self, kk: &mut Tensor2, _phi: f64) -> Result<(), StrError> {
        kk.sym_set(0, 0, self.kx);
        kk.sym_set(1, 1, self.ky);
        if !self.two_dim {
            kk.sym_set(2, 2, self.kz);
        }
        Ok(())
    }
}

// IsotropicLinearConductivityModel ////////////////////////////////////////////////////////////////////////////////////

/// IsotropicLinearConductivityModel
///
/// k = kᵣ·(1 + β·T)
pub struct IsotropicLinearConductivityModel {
    kr: f64,
    beta: f64,
    two_dim: bool,
}

impl IsotropicLinearConductivityModel {
    pub fn new(kr: f64, beta: f64, two_dim: bool) -> Self {
        IsotropicLinearConductivityModel { kr, beta, two_dim }
    }
}

impl ConductivityModel for IsotropicLinearConductivityModel {
    fn tensor(&mut self, kk: &mut Tensor2, tt: f64) -> Result<(), StrError> {
        let val = self.kr * (1.0 * self.beta * tt);
        kk.sym_set(0, 0, val);
        kk.sym_set(1, 1, val);
        if !self.two_dim {
            kk.sym_set(2, 2, val);
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {

    #[test]
    fn allocate_conductivity_model_works() {}
}
