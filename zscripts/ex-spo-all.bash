#!/bin/bash

set -e

echo
echo "============================================="
echo "  Running SPO 751: Pressurized Cylinder"
echo "============================================="
bash "$(dirname "$0")/ex-spo-751_pres_cylin.bash"

echo
echo "============================================="
echo "  Running SPO 752: Pressurized Sphere"
echo "============================================="
bash "$(dirname "$0")/ex-spo-752_pres_sphere.bash"

echo
echo "============================================="
echo "  Running SPO 753: Circular Plate"
echo "============================================="
bash "$(dirname "$0")/ex-spo-753_circ_plate.bash"

echo
echo "============================================="
echo "  Running SPO 754: Footing"
echo "============================================="
bash "$(dirname "$0")/ex-spo-754_footing.bash"

echo
echo "============================================="
echo "  Running SPO 755: Tensile Test"
echo "============================================="
bash "$(dirname "$0")/ex-spo-755_tensile.bash"

echo
echo
echo "✨✨✨✨✨ All SPO examples done! ✨✨✨✨✨"
echo