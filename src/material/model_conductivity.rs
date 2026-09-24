use crate::base::Conductivity;
use crate::StrError;
use crate::D3;
use russell_tensor::Tensor2;

/// Implements conductivity models
///
/// `N` specifies the space dimension and must be [crate::D2] or [crate::D3].
/// It is actually `2*ndim` because it defines the tensor representation.
pub struct ModelConductivity<const N: usize> {
    /// Use Constant model
    cte_enabled: bool,

    /// Use Isotropic model
    iso_enabled: bool,

    /// Isotropic model k = (1 + β T) kᵣ I  (I is the identity tensor)
    iso_kr: f64,

    /// Isotropic model k = (1 + β T) kᵣ I  (I is the identity tensor)
    iso_beta: f64,

    /// x-component of the conductivity tensor
    kx: f64,

    /// y-component of the conductivity tensor
    ky: f64,

    /// z-component of the conductivity tensor
    kz: f64,

    /// Pedroso-Zhang-Ehlers (PZE): λ0 parameter
    pze_lambda_0: f64,

    /// Pedroso-Zhang-Ehlers (PZE): λ1 parameter
    pze_lambda_1: f64,

    /// Pedroso-Zhang-Ehlers (PZE): α parameter
    pze_alpha: f64,

    /// Pedroso-Zhang-Ehlers (PZE): β parameter
    pze_beta: f64,
}

impl<const N: usize> ModelConductivity<N> {
    const VALIDATE_N: () = assert!(N == 4 || N == 6, "N must be 4 or 6 (=2*NDIM)");

    /// Allocates a new instance
    pub fn new(param: &Conductivity) -> Result<Self, StrError> {
        let _ = Self::VALIDATE_N;

        match *param {
            Conductivity::Constant { kx, ky, kz } => Ok(ModelConductivity {
                cte_enabled: true,
                iso_enabled: false,
                iso_kr: 0.0,
                iso_beta: 0.0,
                kx,
                ky,
                kz,
                pze_lambda_0: 0.0,
                pze_lambda_1: 0.0,
                pze_alpha: 0.0,
                pze_beta: 0.0,
            }),
            Conductivity::IsotropicLinear { kr, beta } => Ok(ModelConductivity {
                cte_enabled: false,
                iso_enabled: true,
                iso_kr: kr,
                iso_beta: beta,
                kx: 0.0,
                ky: 0.0,
                kz: 0.0,
                pze_lambda_0: 0.0,
                pze_lambda_1: 0.0,
                pze_alpha: 0.0,
                pze_beta: 0.0,
            }),
            Conductivity::PedrosoZhangEhlers {
                kx,
                ky,
                kz,
                lambda_0,
                lambda_1,
                alpha,
                beta,
            } => Ok(ModelConductivity {
                cte_enabled: false,
                iso_enabled: false,
                iso_kr: 0.0,
                iso_beta: 0.0,
                kx,
                ky,
                kz,
                pze_lambda_0: lambda_0,
                pze_lambda_1: lambda_1,
                pze_alpha: alpha,
                pze_beta: beta,
            }),
        }
    }

    /// Indicates whether or not the model has a symmetric k
    pub fn has_symmetric_k(&self) -> bool {
        true
    }

    /// Indicates whether the conductivity tensor depends on ϕ or not, thus ∂k/∂ϕ is needed
    pub fn has_variable_k(&self) -> bool {
        !self.cte_enabled
    }

    /// Calculates the conductivity tensor
    pub fn calc_k(&self, k: &mut Tensor2<N>, phi: f64) -> Result<(), StrError> {
        k.clear();
        if self.cte_enabled {
            k.sym_set_std(0, 0, self.kx);
            k.sym_set_std(1, 1, self.ky);
            if N == D3 {
                k.sym_set_std(2, 2, self.kz);
            }
        } else if self.iso_enabled {
            // k = (1 + β T) kᵣ I   (I is the identity tensor)
            let val = (1.0 + self.iso_beta * phi) * self.iso_kr;
            k.sym_set_std(0, 0, val);
            k.sym_set_std(1, 1, val);
            if N == D3 {
                k.sym_set_std(2, 2, val);
            }
        } else {
            let _ = self.pze_lambda_0;
            let _ = self.pze_lambda_1;
            let _ = self.pze_alpha;
            let _ = self.pze_beta;
            return Err("TODO: Pedroso-Zhang-Ehlers Conductivity model");
        }
        Ok(())
    }

    /// Calculates the derivative of the conductivity tensor with respect to phi (∂k/∂ϕ)
    pub fn calc_dk_dphi(&self, dk_dphi: &mut Tensor2<N>, _phi: f64) -> Result<(), StrError> {
        dk_dphi.clear();
        if self.cte_enabled {
            // nothing else needed
        } else if self.iso_enabled {
            // k = (1 + β T) kᵣ I   →  ∂k/∂ϕ = ∂k/∂T = β kᵣ I
            let val = self.iso_beta * self.iso_kr;
            dk_dphi.sym_set_std(0, 0, val);
            dk_dphi.sym_set_std(1, 1, val);
            if N == D3 {
                dk_dphi.sym_set_std(2, 2, val);
            }
        } else {
            return Err("TODO: Pedroso-Zhang-Ehlers Conductivity model");
        }
        Ok(())
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#[cfg(test)]
mod tests {
    use super::ModelConductivity;
    use crate::base::Conductivity;
    use crate::D2;
    use russell_lab::{approx_eq, deriv1_central5};
    use russell_tensor::Tensor2;

    #[test]
    fn derivative_works() {
        let param = Conductivity::IsotropicLinear { kr: 20.0, beta: 0.5 };
        let model = ModelConductivity::<D2>::new(&param).unwrap();

        let phi_ini = 100.0;
        let mut dk_dphi_ana = Tensor2::<D2>::new();
        model.calc_dk_dphi(&mut dk_dphi_ana, phi_ini).unwrap();

        struct Args {
            temp: Tensor2<D2>,
        }
        let mut args = Args { temp: Tensor2::new() };

        for i in 0..2 {
            for j in 0..2 {
                let num = deriv1_central5(phi_ini, &mut args, |phi_at, a| {
                    model.calc_k(&mut a.temp, phi_at).unwrap();
                    Ok(a.temp.get_std(i, j))
                })
                .unwrap();
                // println!("k[{},{}] = {:?} → {:?}", i, j, dk_dphi_ana.get(i, j), num);
                approx_eq(dk_dphi_ana.get_std(i, j), num, 1e-10);
            }
        }
    }
}
