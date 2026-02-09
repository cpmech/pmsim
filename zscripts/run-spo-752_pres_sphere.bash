#!/bin/bash

set -e

cargo build --release

cargo run --release --example spo_752_pres_sphere -- -g mumps
cargo run --release --example spo_752_pres_sphere -- -g mumps --lmm
cargo run --release --example spo_752_pres_sphere -- -g mumps --arclength
cargo run --release --example spo_752_pres_sphere -- -g mumps --arclength --lmm

cargo run --release --example spo_752_pres_sphere -- -g mumps --residual  
cargo run --release --example spo_752_pres_sphere -- -g mumps --residual --lmm
cargo run --release --example spo_752_pres_sphere -- -g mumps --residual --arclength
cargo run --release --example spo_752_pres_sphere -- -g mumps --residual --arclength --lmm

cargo run --release --example spo_752_pres_sphere -- -g klu --residual --arclength
cargo run --release --example spo_752_pres_sphere -- -g umfpack --residual --arclength
cargo run --release --example spo_752_pres_sphere -- -g umfpack --residual --arclength --lmm
