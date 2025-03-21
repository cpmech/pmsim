#!/bin/bash

cargo test --test test_spo_752_pres_sphere -- --nocapture
cargo run --bin pmsim_to_paraview -- /tmp/pmsim spo_752_pres_sphere_residual
