#!/bin/bash

set -e

cargo build
cargo test test_spo_755_tensile

~/rust_modules/debug/pmsim_to_paraview /tmp/pmsim/results test_spo_755_tensile
