#!/bin/bash

TESTS=(
    test_spo_751_pres_cylin
    test_spo_752_pres_sphere
    test_spo_753_circ_plate
    test_spo_754_footing
    test_spo_755_tensile
)

for test in "${TESTS[@]}"; do
    cargo test --test "$test" -- --nocapture
done

DIRS=(
    spo_751
    spo_752
    spo_753
    spo_754
    spo_755
)

for dir in "${DIRS[@]}"; do
    cp /tmp/pmsim/"$dir"/*.txt ./data/spo/"$dir"/
done
