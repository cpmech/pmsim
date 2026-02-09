#!/bin/bash

set -e

cargo build --release

cargo run --release --example spo_754_footing -- -g mumps
cargo run --release --example spo_754_footing -- -g mumps --lmm
cargo run --release --example spo_754_footing -- -g mumps --arclength
cargo run --release --example spo_754_footing -- -g mumps --arclength --lmm

cargo run --release --example spo_754_footing -- -g klu --arclength
cargo run --release --example spo_754_footing -- -g umfpack --arclength
cargo run --release --example spo_754_footing -- -g umfpack --arclength --lmm
