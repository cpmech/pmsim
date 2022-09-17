use crate::base::ParamConductivity;
use crate::StrError;
use russell_tensor::Tensor2;

/// Implements conductivity models
pub struct ConductivityModel<'a> {
    /// Material parameters
    param: &'a ParamConductivity,

    /// Indicates a 2D conductivity tensor
    two_dim: bool,

    /// Indicates whether the conductivity tensor depends on ϕ or not
    variable_k: bool,
}

impl<'a> ConductivityModel<'a> {
    /// Allocates a new instance
    pub fn new(param: &'a ParamConductivity, two_dim: bool) -> Self {
        let variable_k = match param {
            ParamConductivity::Constant { .. } => false,
            ParamConductivity::IsotropicLinear { .. } => true,
            ParamConductivity::PedrosoZhangEhlers { .. } => true,
        };
        ConductivityModel {
            param,
            two_dim,
            variable_k,
        }
    }

    /// Calculates the conductivity tensor
    pub fn calc_k(&self, k: &mut Tensor2, phi: f64) -> Result<(), StrError> {
        k.clear();
        match self.param {
            ParamConductivity::Constant { kx, ky, kz } => {
                k.sym_set(0, 0, *kx);
                k.sym_set(1, 1, *ky);
                if !self.two_dim {
                    k.sym_set(2, 2, *kz);
                }
            }
            ParamConductivity::IsotropicLinear { kr, beta } => {
                // k = (1 + β T) kᵣ I   (I is the identity tensor)
                let val = (1.0 + beta * phi) * kr;
                k.sym_set(0, 0, val);
                k.sym_set(1, 1, val);
                if !self.two_dim {
                    k.sym_set(2, 2, val);
                }
            }
            _ => panic!("todo"),
        }
        Ok(())
    }

    /// Indicates whether or not the model has a variable k, thus ∂k/∂ϕ is needed
    pub fn has_variable_k(&self) -> bool {
        self.variable_k
    }

    /// Calculates the derivative of the conductivity tensor with respect to phi (∂k/∂ϕ)
    pub fn calc_dk_dphi(&self, dk_dphi: &mut Tensor2, _phi: f64) -> Result<(), StrError> {
        dk_dphi.clear();
        match self.param {
            ParamConductivity::Constant { .. } => (),
            ParamConductivity::IsotropicLinear { kr, beta } => {
                // k = (1 + β T) kᵣ I   →  ∂k/∂ϕ = ∂k/∂T = β kᵣ I
                let val = beta * kr;
                dk_dphi.sym_set(0, 0, val);
                dk_dphi.sym_set(1, 1, val);
                if !self.two_dim {
                    dk_dphi.sym_set(2, 2, val);
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
    use super::ConductivityModel;
    use crate::base::ParamConductivity;
    use russell_chk::{approx_eq, deriv_central5};
    use russell_tensor::Tensor2;

    #[test]
    fn derivative_works() {
        let param = ParamConductivity::IsotropicLinear { kr: 20.0, beta: 0.5 };
        let model = ConductivityModel::new(&param, true);

        let phi_ini = 100.0;
        let mut dk_dphi_ana = Tensor2::new(true, true);
        model.calc_dk_dphi(&mut dk_dphi_ana, phi_ini).unwrap();

        struct Args {
            temp: Tensor2,
        }
        let mut args = Args {
            temp: Tensor2::new(true, true),
        };

        for i in 0..2 {
            for j in 0..2 {
                let num = deriv_central5(phi_ini, &mut args, |phi_at, a| {
                    model.calc_k(&mut a.temp, phi_at).unwrap();
                    a.temp.get(i, j)
                });
                // println!("k[{},{}] = {:?} → {:?}", i, j, dk_dphi_ana.get(i, j), num);
                approx_eq(dk_dphi_ana.get(i, j), num, 1e-10);
            }
        }
    }
}
