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
run_example spo_751_pres_cylin -- -g mumps
run_example spo_751_pres_cylin -- -g mumps --lmm
run_example spo_751_pres_cylin -- -g mumps --arclength
run_example spo_751_pres_cylin -- -g mumps --arclength --lmm
run_example spo_751_pres_cylin -- -g mumps --arclength --bordering
run_example spo_751_pres_cylin -- -g mumps --arclength --lmm --bordering

# mumps: residual
run_example spo_751_pres_cylin -- -g mumps --residual  
run_example spo_751_pres_cylin -- -g mumps --residual --lmm
run_example spo_751_pres_cylin -- -g mumps --residual --arclength
run_example spo_751_pres_cylin -- -g mumps --residual --arclength --lmm
run_example spo_751_pres_cylin -- -g mumps --residual --arclength --bordering
run_example spo_751_pres_cylin -- -g mumps --residual --arclength --lmm --bordering

# klu: residual
run_example spo_751_pres_cylin -- -g klu --residual --arclength

# umfpack: residual
run_example spo_751_pres_cylin -- -g umfpack --residual --arclength
run_example spo_751_pres_cylin -- -g umfpack --residual --arclength --lmm
run_example spo_751_pres_cylin -- -g umfpack --residual --arclength --bordering
run_example spo_751_pres_cylin -- -g umfpack --residual --arclength --lmm --bordering

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo