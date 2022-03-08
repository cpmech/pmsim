use crate::{ParamConductivity, StrError};
use russell_tensor::Tensor2;

/// Implements a model for liquid conductivity within porous media
///
/// # Reference
///
/// * Pedroso DM, Zhang Y, Ehlers W (2017) Solution of liquid-gas-solid coupled
///   equations for porous media considering dynamics and hysteretic behavior,
///   ASCE Journal of Engineering Mechanics, 143:6(04017021),
///   <https://dx.doi.org/10.1061/(ASCE)EM.1943-7889.0001208>
pub struct ModelConductivity {
    // for constant, linear or pedroso-zhang-ehlers models
    kk_sat: Tensor2,

    // for constant model
    is_constant: bool,

    // for linear model
    is_linear: bool,
    lambda: f64,

    // for pedroso-zhang-ehlers model
    lambda_0: f64,
    lambda_1: f64,
    alpha: f64,
    beta: f64,
    b: f64,
    c1: f64,
    c2: f64,
    c3: f64,
}

impl ModelConductivity {
    /// Allocates a new instance
    pub fn new(param: &ParamConductivity, two_dim: bool) -> Result<Self, StrError> {
        let mut model = ModelConductivity {
            kk_sat: Tensor2::new(true, two_dim),
            is_constant: false,
            is_linear: false,
            lambda: 0.0,
            lambda_0: 0.0,
            lambda_1: 0.0,
            alpha: 0.0,
            beta: 0.0,
            b: 0.0,
            c1: 0.0,
            c2: 0.0,
            c3: 0.0,
        };
        let (kx, ky, kz) = match param {
            &ParamConductivity::Constant { kx, ky, kz } => {
                model.is_constant = true;
                (kx, ky, kz)
            }
            &ParamConductivity::Linear { kx, ky, kz, lambda } => {
                if lambda < 1.0 {
                    return Err("lambda must be greater than or equal to 1.0 for linear model");
                }
                model.is_linear = true;
                model.lambda = lambda;
                (kx, ky, kz)
            }
            &ParamConductivity::PedrosoZhangEhlers {
                kx,
                ky,
                kz,
                lambda_0,
                lambda_1,
                alpha,
                beta,
            } => {
                let b = if lambda_1 < lambda_0 { 1.0 } else { -1.0 };
                model.lambda_0 = lambda_0;
                model.lambda_1 = lambda_1;
                model.alpha = alpha;
                model.beta = beta;
                model.b = b;
                model.c1 = beta * b * (lambda_1 - lambda_0);
                model.c2 = f64::exp(beta * b * alpha);
                model.c3 = f64::exp(beta * b * (1.0 - lambda_0)) - model.c2 * f64::exp(model.c1);
                (kx, ky, kz)
            }
        };
        if kx < 0.0 {
            return Err("kx must be greater than zero");
        }
        if ky < 0.0 {
            return Err("ky must be greater than zero");
        }
        model.kk_sat.sym_set(0, 0, kx);
        model.kk_sat.sym_set(1, 1, ky);
        if !two_dim {
            if kz < 0.0 {
                return Err("kz must be greater than zero");
            }
            model.kk_sat.sym_set(2, 2, kz);
        }
        Ok(model)
    }

    /// Returns an access to the saturated conductivity tensor
    pub fn conductivity_tensor(&self) -> &Tensor2 {
        &self.kk_sat
    }

    /// Returns the relative_conductivity for given saturation
    pub fn relative_conductivity(&self, saturation: f64) -> Result<f64, StrError> {
        if saturation < 0.0 || saturation > 1.0 {
            return Err("saturation must be in 0 ≤ saturation ≤ 1 to calculate conductivity");
        }
        if self.is_constant {
            return Ok(1.0);
        }
        if self.is_linear {
            if self.lambda * saturation > 1.0 {
                return Ok(1.0);
            } else {
                return Ok(self.lambda * saturation);
            };
        }
        let kr = self.lambda_0 * saturation
            + f64::ln(self.c3 + self.c2 * f64::exp(self.c1 * saturation)) / (self.beta * self.b);
        return Ok(kr);
    }
}
