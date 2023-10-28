#!/bin/bash

set -e

for ex in examples/*.rs; do
    key=`basename -s ".rs" $ex`
    if [ "$key" = "pressurized_cylinder_plot" ] || 
       [ "$key" = "pressurized_cylinder_table" ] || 
       [ "$key" = "pressurized_cylinder2d_elast" ] || 
       [ "$key" = "pressurized_cylinder3d_elast" ]; then
        echo
        echo "skip example $key"
        echo
    else
        cargo run --example $key
    fi
done
