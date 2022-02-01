#![allow(dead_code, unused_mut, unused_variables)]

use crate::ElementConfig;
use std::collections::HashMap;

pub struct Configuration {
    element_config: HashMap<usize, ElementConfig>,
}

impl Configuration {
    pub fn new() -> Self {
        Configuration {
            element_config: HashMap::new(),
        }
    }
}
