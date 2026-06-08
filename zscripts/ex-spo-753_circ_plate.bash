#!/bin/bash

set -e

# define features
FEAT="--features intel_mkl,local_sparse"

# build the example
cargo build --release $FEAT

# run the examples with options
RUN="cargo run --release $FEAT --example"

# mumps
$RUN spo_753_circ_plate -- -g mumps
#$RUN spo_753_circ_plate -- -g mumps --lmm
$RUN spo_753_circ_plate -- -g mumps --arclength
$RUN spo_753_circ_plate -- -g mumps --arclength --lmm
$RUN spo_753_circ_plate -- -g mumps --arclength --bordering
$RUN spo_753_circ_plate -- -g mumps --arclength --lmm --bordering

# klu
$RUN spo_753_circ_plate -- -g klu --arclength

# umfpack
$RUN spo_753_circ_plate -- -g umfpack
$RUN spo_753_circ_plate -- -g umfpack --arclength
$RUN spo_753_circ_plate -- -g umfpack --arclength --bordering
$RUN spo_753_circ_plate -- -g umfpack --arclength --lmm
$RUN spo_753_circ_plate -- -g umfpack --arclength --lmm --bordering

# comparison
$RUN spo_753_circ_plate_compare -- -g umfpack --arclength --bord-vs-full
$RUN spo_753_circ_plate_compare -- -g umfpack --arclength --lmm --bord-vs-full
$RUN spo_753_circ_plate_compare -- -g umfpack --arclength --lmm-vs-sps

echo
echo
echo "✨✨✨✨✨ All done! ✨✨✨✨✨"
echo