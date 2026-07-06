#!/bin/bash

TESTS="
test_heat_arpaci_nonlinear_1d
test_von_mises_single_element_2d
test_von_mises_2x2_elements_2d
test_spo_751_pres_cylin
"

for test in $TESTS; do
    echo
    echo
    echo
    cargo test --package pmsim --test $test -- $test --exact --show-output
done
