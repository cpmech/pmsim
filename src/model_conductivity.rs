use crate::{ParamConductivity, StrError};

pub struct ModelConductivity {
    // for constant, linear or pedroso-zhang-ehlers models
    kx: f64,
    ky: f64,
    kz: f64,
    is_constant: bool,

    // for linear model
    lambda: f64,
    is_linear: bool,

    // for pedroso-zhang-ehlers model
    lambda_0: f64,
    lambda_1: f64,
    alpha: f64,
    beta: f64,
}

impl ModelConductivity {
    pub fn new(params: &ParamConductivity) -> Result<Self, StrError> {
        let mut model = ModelConductivity {
            // for constant, linear or pedroso-zhang-ehlers models
            kx: 0.0,
            ky: 0.0,
            kz: 0.0,
            is_constant: false,
            // for linear model
            lambda: 0.0,
            is_linear: false,
            // for pedroso-zhang-ehlers model
            lambda_0: 0.0,
            lambda_1: 0.0,
            alpha: 0.0,
            beta: 0.0,
        };

        match params {
            &ParamConductivity::Constant { kx, ky, kz } => {
                model.kx = kx;
                model.ky = ky;
                model.kz = kz;
                model.is_constant = true;
            }
            &ParamConductivity::Linear { kx, ky, kz, lambda } => {
                model.kx = kx;
                model.ky = ky;
                model.kz = kz;
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
                model.kx = kx;
                model.ky = ky;
                model.kz = kz;
                model.lambda_0 = lambda_0;
                model.lambda_1 = lambda_1;
                model.alpha = alpha;
                model.beta = beta;
            }
        };

        // todo: check parameters

        Ok(model)
    }
}
