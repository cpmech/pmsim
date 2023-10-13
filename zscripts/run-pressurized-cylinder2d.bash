#!/bin/bash

GENIES="\
    mumps \
    umfpack \
    inteldss \
"

KINDS="\
    tri3 \
    tri6 \
    tri10 \
    tri15 \
    qua4 \
    qua8 \
    qua9 \
    qua12 \
    qua16 \
    qua17 \
"

for genie in $GENIES; do
    echo
    for kind in $KINDS; do
        cargo run --release --example pressurized_cylinder2d_elast -- $genie $kind
    done
    cargo run --release --example pressurized_cylinder2d_elast_results -- $genie
done
