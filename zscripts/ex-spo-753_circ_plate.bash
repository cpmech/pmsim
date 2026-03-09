#!/bin/bash

set -e

cargo build --release

# cargo run --release --example spo_753_circ_plate -- -g mumps
# cargo run --release --example spo_753_circ_plate -- -g mumps --lmm
# cargo run --release --example spo_753_circ_plate -- -g mumps --arclength
# cargo run --release --example spo_753_circ_plate -- -g mumps --arclength --lmm
# cargo run --release --example spo_753_circ_plate -- -g mumps --arclength --bordering
# cargo run --release --example spo_753_circ_plate -- -g mumps --arclength --lmm --bordering

#cargo run --release --example spo_753_circ_plate -- -g klu --arclength
#cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength
#cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength --bordering

# cargo run --release --example spo_753_circ_plate -- -g umfpack
cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength
# cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength --bordering
# cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength --lmm
# cargo run --release --example spo_753_circ_plate -- -g umfpack --arclength --lmm --bordering

# cargo run --release --example spo_753_circ_plate_compare -- -g umfpack --arclength --bord-vs-full
# cargo run --release --example spo_753_circ_plate_compare -- -g umfpack --arclength --lmm --bord-vs-full
# cargo run --release --example spo_753_circ_plate_compare -- -g umfpack --arclength --lmm-vs-sps

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo