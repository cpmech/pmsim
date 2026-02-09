#!/bin/bash

set -e

cargo build --release

cargo run --release --example spo_755_tensile -- -g mumps
cargo run --release --example spo_755_tensile -- -g mumps --lmm
cargo run --release --example spo_755_tensile -- -g mumps --arclength
cargo run --release --example spo_755_tensile -- -g mumps --arclength --lmm

cargo run --release --example spo_755_tensile -- -g klu --arclength
cargo run --release --example spo_755_tensile -- -g umfpack --arclength
cargo run --release --example spo_755_tensile -- -g umfpack --arclength --lmm

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo