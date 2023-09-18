#!/bin/bash

cargo run --release --example pressurized_cylinder2d_elast -- tri3
cargo run --release --example pressurized_cylinder2d_elast -- tri6
cargo run --release --example pressurized_cylinder2d_elast -- tri10
cargo run --release --example pressurized_cylinder2d_elast -- tri15
cargo run --release --example pressurized_cylinder2d_elast -- qua4
cargo run --release --example pressurized_cylinder2d_elast -- qua8
cargo run --release --example pressurized_cylinder2d_elast -- qua9
cargo run --release --example pressurized_cylinder2d_elast -- qua12
cargo run --release --example pressurized_cylinder2d_elast -- qua16
cargo run --release --example pressurized_cylinder2d_elast -- qua17

cargo run --release --example pressurized_cylinder2d_elast_results
