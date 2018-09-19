#!/usr/bin/env bash

module load python

export CHARGE_ACCOUNT=cli106ice

# -----------------
# GENERATE BASELINE
# -----------------
./create_test MVK_PL.ne4_ne4.FC5AV1C-04P2 -p cli106 --baseline-root "${PROJWORK}/cli106/kennedy/baselines" -g -o --walltime 02:00:00

# -----------------
# COMPARE BASELINE
# -----------------
# NOTE: -C is required!
# -----------------
# ./create_test MVK_PL.ne4_ne4.FC5AV1C-04P2 -p cli106 --baseline-root "${PROJWORK}/cli106/kennedy/baselines" -c --walltime 02:00:00
