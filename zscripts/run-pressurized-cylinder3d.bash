#!/bin/bash

GENIES="\
    umfpack \
"

KINDS="\
    hex8 \
    hex20 \
    tet4 \
    tet10
"

# KINDS="tet10"

for genie in $GENIES; do
    echo
    for kind in $KINDS; do
        cargo run --release --example pressurized_cylinder3d_elast -- $genie $kind
    done
    cargo run --release --example plot_convergence_results -- pressurized_cylinder3d_elast $genie
done
