use russell_tensor::Tensor2;

use crate::{ParamConductivity, StrError};

pub struct ModelConductivity {
    // for constant, linear or pedroso-zhang-ehlers models
    kk_sat: Tensor2,
    is_constant: bool,

    // for linear model
    lambda: f64,
    is_linear: bool,

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
    pub fn new(params: &ParamConductivity, two_dim: bool) -> Result<Self, StrError> {
        let mut model = ModelConductivity {
            // for constant, linear or pedroso-zhang-ehlers models
            kk_sat: Tensor2::new(true, two_dim),
            is_constant: false,
            // for linear model
            lambda: 0.0,
            is_linear: false,
            // for pedroso-zhang-ehlers model
            lambda_0: 0.0,
            lambda_1: 0.0,
            alpha: 0.0,
            beta: 0.0,
            b: 0.0,
            c1: 0.0,
            c2: 0.0,
            c3: 0.0,
        };

        match params {
            &ParamConductivity::Constant { kx, ky, kz } => {
                model.kk_sat.sym_set(0, 0, kx);
                model.kk_sat.sym_set(1, 1, ky);
                if !two_dim {
                    model.kk_sat.sym_set(2, 2, kz);
                }
                model.is_constant = true;
            }
            &ParamConductivity::Linear { kx, ky, kz, lambda } => {
                model.kk_sat.sym_set(0, 0, kx);
                model.kk_sat.sym_set(1, 1, ky);
                if !two_dim {
                    model.kk_sat.sym_set(2, 2, kz);
                }
                model.lambda = lambda;
                model.is_linear = true;
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
                model.kk_sat.sym_set(0, 0, kx);
                model.kk_sat.sym_set(1, 1, ky);
                if !two_dim {
                    model.kk_sat.sym_set(2, 2, kz);
                }
                let b = if lambda_1 < lambda_0 { 1.0 } else { -1.0 };
                model.lambda_0 = lambda_0;
                model.lambda_1 = lambda_1;
                model.alpha = alpha;
                model.beta = beta;
                model.b = b;
                model.c1 = beta * b * (lambda_1 - lambda_0);
                model.c2 = f64::exp(beta * b * alpha);
                model.c3 = f64::exp(beta * b * (1.0 - lambda_0)) - model.c2 * f64::exp(model.c1);
            }
        };

        // todo: check parameters

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
