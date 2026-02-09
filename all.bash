#!/bin/bash

cargo test

bash ./zscripts/run-spo-751_pres_cylin.bash
bash ./zscripts/run-spo-752_pres_sphere.bash
