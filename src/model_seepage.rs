#![allow(dead_code, unused_mut, unused_variables, unused_imports)]

use crate::{ParamSeepageLiq, StrError};

pub struct ModelSeepageLiq {
    // todo
}

impl ModelSeepageLiq {
    pub fn new(params: &ParamSeepageLiq) -> Result<Self, StrError> {
        let model = ModelSeepageLiq {};
        Ok(model)
    }
}
