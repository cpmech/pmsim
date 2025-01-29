#!/bin/bash

set -e

cargo build
cargo test test_spo_751_press_cylin
~/rust_modules/debug/pmsim_to_paraview /tmp/pmsim/results test_spo_751_pres_cylin
