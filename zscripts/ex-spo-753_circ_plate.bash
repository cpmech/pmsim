#!/bin/bash

set -e

cargo build --release

cargo run --release --example spo_753_circ_plate -- -g mumps
cargo run --release --example spo_753_circ_plate -- -g mumps --lmm
cargo run --release --example spo_753_circ_plate -- -g mumps --arclength
cargo run --release --example spo_753_circ_plate -- -g mumps --arclength --lmm

cargo run --release --example spo_753_circ_plate -- -g klu --arclength
cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength
cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength --lmm

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo