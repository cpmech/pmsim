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

TOL=1e-9
run_example spo_753_circ_plate -- -g cudss --arclength --tol $TOL
run_example spo_753_circ_plate -- -g mumps --arclength --tol $TOL
run_example spo_753_circ_plate -- -g umfpack --arclength --tol $TOL
run_example spo_753_circ_plate_compare_genies -- --arclength --tol $TOL

TOL=1e-10
run_example spo_753_circ_plate -- -g cudss --arclength --tol $TOL
run_example spo_753_circ_plate -- -g mumps --arclength --tol $TOL
run_example spo_753_circ_plate -- -g umfpack --arclength --tol $TOL
run_example spo_753_circ_plate_compare_genies -- --arclength --tol $TOL

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo