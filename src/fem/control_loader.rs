use super::FemState;
use crate::base::Config;
use crate::StrError;

/// Implements the loading control using the lambda parameter
pub(crate) struct ControlLoader<'a> {
    config: &'a Config<'a>,
    last: bool,
}

impl<'a> ControlLoader<'a> {
    /// Allocates a new instance
    pub fn new(config: &'a Config) -> Self {
        ControlLoader { config, last: false }
    }

    pub fn initialize(&mut self, state: &mut FemState) {
        state.lambda = 0.0;
        self.last = false;
    }

    /// Returns whether the last (time) loading increment (lambda) has been reached
    pub fn last(&self) -> bool {
        self.last
    }

    pub fn next(&mut self, state: &mut FemState) -> Result<(), StrError> {
        // set first Δλ
        if state.step == 0 {
            state.ddl = self.config.ddl;
        }

        // check for Δλ too small
        if state.ddl < self.config.ddl_min {
            return Err("Δλ is smaller than the allowed minimum");
        }

        // check for final loading increment
        if state.lambda + state.ddl >= 1.0 {
            if !self.config.constant_ddl && state.lambda + state.ddl != 1.0 {
                // only truncates if λ+Δλ is not exactly equal to 1.0
                state.ddl = f64::max(self.config.ddl_min, 1.0 - state.lambda);
            }
            self.last = true;
        }

        // update λ
        state.lambda += state.ddl;

        // trial displacement u
        // the trial displacement is the previous displacement → do nothing

        // displacement increment Δu
        state.ddu.fill(0.0);
        Ok(())
    }
}
