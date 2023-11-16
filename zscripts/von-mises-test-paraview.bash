#!/bin/bash

DIR="/tmp/pmsim/results"
STEM="test_von_mises_single_element_2d"

cargo test --test $STEM
cargo run --bin pmsim_to_paraview -- $DIR $STEM
