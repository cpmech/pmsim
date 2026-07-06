#!/bin/bash

set -e

# define features
FEAT="--all-features"

# build the example
cargo build --release $FEAT

# run the examples with options
run_example() {
    echo
    echo "cargo run --release $FEAT --example $*"
    echo
    cargo run --release $FEAT --example "$@"
}

# mumps
run_example spo_753_circ_plate -- -g mumps
run_example spo_753_circ_plate -- -g mumps --lmm
run_example spo_753_circ_plate -- -g mumps --arclength
run_example spo_753_circ_plate -- -g mumps --arclength --lmm
run_example spo_753_circ_plate -- -g mumps --arclength --bordering
run_example spo_753_circ_plate -- -g mumps --arclength --lmm --bordering

# cudss
run_example spo_753_circ_plate -- -g cudss --arclength

# umfpack
run_example spo_753_circ_plate -- -g umfpack
run_example spo_753_circ_plate -- -g umfpack --lmm
run_example spo_753_circ_plate -- -g umfpack --arclength
run_example spo_753_circ_plate -- -g umfpack --arclength --bordering
run_example spo_753_circ_plate -- -g umfpack --arclength --lmm
run_example spo_753_circ_plate -- -g umfpack --arclength --lmm --bordering

# comparison
run_example spo_753_circ_plate_compare -- -g umfpack --arclength --bord-vs-full
run_example spo_753_circ_plate_compare -- -g umfpack --arclength --lmm --bord-vs-full
run_example spo_753_circ_plate_compare -- -g umfpack --arclength --lmm-vs-sps

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo