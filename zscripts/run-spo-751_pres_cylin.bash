#!/bin/bash

set -e

cargo build --release

# cargo run --release --example spo_751_pres_cylin -- -g mumps
# cargo run --release --example spo_751_pres_cylin -- -g mumps --lmm
# cargo run --release --example spo_751_pres_cylin -- -g mumps --arclength
# cargo run --release --example spo_751_pres_cylin -- -g mumps --arclength --lmm

# cargo run --release --example spo_751_pres_cylin -- -g mumps --residual  
# cargo run --release --example spo_751_pres_cylin -- -g mumps --residual --lmm
# cargo run --release --example spo_751_pres_cylin -- -g mumps --residual --arclength
cargo run --release --example spo_751_pres_cylin -- -g mumps --residual --arclength --lmm