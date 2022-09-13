use crate::base::ParamConductivity;
use crate::StrError;
use russell_tensor::Tensor2;

pub struct ConductivityModel<'a> {
    param: &'a ParamConductivity,
    two_dim: bool,
}

impl<'a> ConductivityModel<'a> {
    pub fn new(param: &'a ParamConductivity, two_dim: bool) -> Self {
        ConductivityModel { param, two_dim }
    }

    pub fn tensor(&mut self, kk: &mut Tensor2, tt: f64) -> Result<(), StrError> {
        match self.param {
            ParamConductivity::Constant { kx, ky, kz } => {
                kk.sym_set(0, 0, *kx);
                kk.sym_set(1, 1, *ky);
                if !self.two_dim {
                    kk.sym_set(2, 2, *kz);
                }
            }
            ParamConductivity::IsotropicLinear { kr, beta } => {
                // k = kᵣ·(1 + β·T)
                let val = kr * (1.0 + beta * tt);
                kk.sym_set(0, 0, val);
                kk.sym_set(1, 1, val);
                if !self.two_dim {
                    kk.sym_set(2, 2, val);
                }
            }
            _ => panic!("todo"),
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
