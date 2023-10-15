#!/bin/bash

set -e

for ex in examples/*.rs; do
    key=`basename -s ".rs" $ex`
    if [ "$key" = "pressurized_cylinder2d_elast" ] || [ "$key" = "pressurized_cylinder2d_elast_results" ]; then
        echo
        echo "skip example $key"
        echo
    else
        cargo run --example $key
    fi
done
