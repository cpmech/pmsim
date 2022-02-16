use crate::{ModelBrooksCorey, ModelPedrosoWilliams, ModelVanGenuchten, ParamLiquidRetention, StrError};
use russell_lab::Vector;

/// Defines a trait for models for liquid retention in porous media
pub trait ModelLiquidRetentionTrait {
    /// Returns the saturation limits (sl_min,sl_max)
    fn saturation_limits(&self) -> (f64, f64);

    /// Calculates Cc(pc,sl) = ∂sl/∂pc
    fn calc_cc(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError>;

    /// Calculates J = dCc/dsl
    fn calc_dcc_dsl(&self, pc: f64, sl: f64, wetting: bool) -> Result<f64, StrError>;
}

/// Generalizes a model for liquid retention in porous media
pub struct ModelLiquidRetention {
    pub model: Box<dyn ModelLiquidRetentionTrait>,
    update_nit_max: usize, // max number of iterations for the update_saturation function
    update_tolerance: f64, // tolerance for the update_saturation function
}

impl ModelLiquidRetention {
    /// Returns a new instance of ModelLiquidRetention
    pub fn new(params: &ParamLiquidRetention) -> Result<Self, StrError> {
        let model: Box<dyn ModelLiquidRetentionTrait> = match params {
            &ParamLiquidRetention::BrooksCorey {
                lambda,
                pc_ae,
                sl_min,
                sl_max,
            } => Box::new(ModelBrooksCorey::new(lambda, pc_ae, sl_min, sl_max)?),
            &ParamLiquidRetention::VanGenuchten {
                alpha,
                m,
                n,
                sl_min,
                sl_max,
                pc_min,
            } => Box::new(ModelVanGenuchten::new(alpha, m, n, sl_min, sl_max, pc_min)?),
            &ParamLiquidRetention::PedrosoWilliams {
                with_hysteresis,
                lambda_d,
                lambda_w,
                beta_d,
                beta_w,
                beta_1,
                beta_2,
                x_rd,
                x_rw,
                y_0,
                y_r,
            } => Box::new(ModelPedrosoWilliams::new(
                with_hysteresis,
                lambda_d,
                lambda_w,
                beta_d,
                beta_w,
                beta_1,
                beta_2,
                x_rd,
                x_rw,
                y_0,
                y_r,
            )?),
        };
        Ok(ModelLiquidRetention {
            model,
            update_nit_max: 12,
            update_tolerance: 1.0e-6,
        })
    }

    /// Returns the maximum liquid saturation
    pub fn get_sl_max(&self) -> f64 {
        self.model.saturation_limits().1
    }

    /// Returns the updated saturation sl_new for given (pc,sl) and Δpc = pc_new - pc
    ///
    /// See Algorithm 2 in the following reference:
    ///
    /// * Pedroso DM (2015) A consistent u-p formulation for porous media with hysteresis,
    ///   Int. J. for Numerical Methods in Engineering, 101:606-634, DOI: 10.1002/nme.4808
    pub fn update_saturation(&self, pc: f64, sl: f64, delta_pc: f64) -> Result<f64, StrError> {
        // constants
        let (sl_min, sl_max) = self.model.saturation_limits();
        // wetting flag
        let wetting = delta_pc < 0.0;
        // forward Euler update
        let f_a = self.model.calc_cc(pc, sl, wetting)?;
        let sl_fe = f64::min(f64::max(sl + delta_pc * f_a, sl_min), sl_max);
        // modified Euler update
        let pc_new = pc + delta_pc;
        let f_b = self.model.calc_cc(pc_new, sl_fe, wetting)?;
        let mut sl_new = f64::min(f64::max(sl + 0.5 * delta_pc * (f_a + f_b), sl_min), sl_max);
        let mut converged = false;
        for _ in 0..self.update_nit_max {
            let cc = self.model.calc_cc(pc_new, sl_new, wetting)?;
            let r = sl_new - sl - delta_pc * cc;
            if f64::abs(r) <= self.update_tolerance {
                converged = true;
                break;
            }
            let dcc_dsl = self.model.calc_dcc_dsl(pc_new, sl_new, wetting)?;
            let dsl = -r / (1.0 - delta_pc * dcc_dsl);
            sl_new += dsl;
        }
        if !converged {
            return Err("Newton-Raphson method failed to converge for the (pc,sl) update");
        }
        Ok(sl_new)
    }

    /// Generates curve data for plotting or other needs
    pub fn generate_curve_data(
        &self,
        pc_start: f64,
        pc_stop: f64,
        sl_start: f64,
        npoint: usize,
    ) -> Result<(Vector, Vector), StrError> {
        if npoint < 2 {
            return Err("at least 2 points are required");
        }
        let all_pc = Vector::linspace(pc_start, pc_stop, npoint)?;
        let mut all_sl = Vector::new(npoint);
        all_sl[0] = sl_start;
        for i in 0..npoint - 1 {
            all_sl[i + 1] = self.update_saturation(all_pc[i], all_sl[i], all_pc[i + 1] - all_pc[i])?;
        }
        Ok((all_pc, all_sl))
    }
}
