#!/bin/bash

cargo test --test test_spo_751_pres_cylin -- --nocapture > data/logs/test_spo_751_pres_cylin.txt
cargo test --test test_stepsize_adaptation_0 -- --nocapture > data/logs/test_stepsize_adaptation_0.txt
cargo test --test test_stepsize_adaptation_1 -- --nocapture > data/logs/test_stepsize_adaptation_1.txt
