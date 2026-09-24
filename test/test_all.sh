#!/usr/bin/env bash
#
# Quick end-to-end test of scripts/run_all.sh with few replications, before a full (100K) run.
# From the main project directory:
#
#   test/test_all.sh [run_all.sh options]
#
# Defaults: all NBA seasons, math_elim_mode = -2 (needs Gurobi), 20 replications (50 for the
# sensitivity run), 4 parallel jobs (which also tests splitting the simulation and aggregating it),
# plots on, results in results/test_run. Any option is passed on to scripts/run_all.sh and overrides
# these defaults, e.g.,
#   test/test_all.sh -m 0 -N                  # without Gurobi and without plots
#   test/test_all.sh -j 1                     # without splitting the simulation
#   test/test_all.sh -n 100 -o results/test2  # more replications, another directory
# Afterwards, test/check_results.jl checks all outputs; the test fails if any check fails.
# It takes about 10-20 minutes (most of it Julia startup and the noisy-rankings experiment).

set -euo pipefail
cd "$(dirname "$0")/.."

OUTDIR=results/test_run
# Find -o in the arguments to know where the results go
args=("$@")
for ((i = 0; i < ${#args[@]}; i++)); do
  if [[ ${args[$i]} == -o ]]; then OUTDIR=${args[$((i+1))]}; fi
done

rm -rf "$OUTDIR"
start=$(date +%s)
scripts/run_all.sh -o "$OUTDIR" -n 20 -S 50 -s all -j 4 "$@"
echo "run_all.sh took $(( $(date +%s) - start )) seconds"
${JULIA:-julia} --project=. test/check_results.jl "$OUTDIR"
