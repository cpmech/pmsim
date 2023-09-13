#!/bin/bash

for ex in examples/*.rs; do
    key=`basename -s ".rs" $ex`
    cargo run --example $key
done
