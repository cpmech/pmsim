#!/bin/bash

set -e

cargo build --release

cargo run --release --example spo_754_footing
