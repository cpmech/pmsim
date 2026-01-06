#!/bin/bash

TESTS=(
    test_spo_751_pres_cylin
    test_spo_754_footing
    test_spo_755_tensile
    test_spo_von_mises_2x2_elements
    test_spo_von_mises_single_element
    test_von_mises_single_element_2d
)

mkdir -p /tmp/pmsim
for test in "${TESTS[@]}"; do
    cargo test --test "$test" -- --nocapture > "/tmp/pmsim/log-${test}.txt"
done
