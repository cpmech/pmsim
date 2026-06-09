#!/bin/bash

set -e

# define features
FEAT="--features intel_mkl,local_sparse"

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
run_example spo_754_footing -- -g mumps
run_example spo_754_footing -- -g mumps --lmm
run_example spo_754_footing -- -g mumps --arclength
run_example spo_754_footing -- -g mumps --arclength --lmm
run_example spo_754_footing -- -g mumps --bordering
run_example spo_754_footing -- -g mumps --lmm --bordering
run_example spo_754_footing -- -g mumps --arclength --bordering
run_example spo_754_footing -- -g mumps --arclength --lmm --bordering

# klu
run_example spo_754_footing -- -g klu --arclength
run_example spo_754_footing -- -g klu --arclength --bordering

# umfpack
run_example spo_754_footing -- -g umfpack --arclength
run_example spo_754_footing -- -g umfpack --arclength --lmm
run_example spo_754_footing -- -g umfpack --arclength --bordering
run_example spo_754_footing -- -g umfpack --arclength --lmm --bordering

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo