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
run_example spo_755_tensile -- -g mumps
run_example spo_755_tensile -- -g mumps --lmm
run_example spo_755_tensile -- -g mumps --arclength
run_example spo_755_tensile -- -g mumps --arclength --lmm

# cudss
run_example spo_755_tensile -- -g cudss --arclength

# umfpack
run_example spo_755_tensile -- -g umfpack --arclength
run_example spo_755_tensile -- -g umfpack --arclength --lmm

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo