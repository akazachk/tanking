#!/usr/bin/env bash
#
# Recreate all experiments and plots, from the main project directory:
#
#   scripts/run_all.sh [options]
#
# Options
#   -o DIR     results directory (default: results/<yyyy-mm-dd>, the date when the run starts)
#   -n N       replications for model validation, the simulation, and noisy rankings (default: 100000)
#   -S N       replications for the sensitivity run (default: 10000)
#   -s SEASONS NBA seasons for model validation, NBA parsing, and Bradley-Terry: all (default) or 2004-2019
#   -g GAMMA   gamma for the simulation and sensitivity run: "auto" (default), the value chosen by model validation
#              (minimax rule; needs the validate experiment, or DIR/gamma.txt from an earlier run), or a number
#   -m MODE    math_elim_mode for the simulation (default: -2, which needs Gurobi)
#   -j JOBS    number of parallel jobs (default: 1); with JOBS > 1, each model in model validation
#              (Bradley-Terry and each gamma) and each of the 31 steps of the simulation is run in its
#              own process, and the results are then combined
#   -t THREADS threads for the sensitivity run (default: the number of jobs, or 4 if JOBS = 1)
#   -e LIST    comma-separated experiments to run, from
#              validate,simulate,parse,noisy,sensitivity,bt (default: validate,simulate,parse,noisy,sensitivity)
#   -N         do not create plots (plots need PyPlot and LaTeX)
#
# Environment variables
#   JULIA      julia command (default: julia); must be Julia 1.6 (the pinned packages do not work with later
#              versions), e.g., JULIA="julia +1.6" with juliaup, or JULIA="julia --sysimage=build/JuliaTanking.so"
#   TANKING_GUROBI_THREADS  threads per Gurobi MIP (default: 1 when JOBS > 1, so parallel jobs do not
#              oversubscribe the cores; otherwise Gurobi's default, all cores)
#
# To use only some cores (e.g., the performance cores), run the script under taskset; all processes it
# starts inherit the CPU affinity:
#   taskset -c 0,2,4,6,8,10,12,14 scripts/run_all.sh -j 8
#
# Output: results in DIR (plots in DIR/pdf and DIR/png), sensitivity results in DIR/sensitivity,
# logs of every step in DIR/logs, and a summary of the settings in DIR/settings.txt.
#
# Examples
#   scripts/run_all.sh -o results/final -j 16                 # everything, all seasons, 100K replications
#   test/test_all.sh                                          # quick test with few replications

set -euo pipefail

cd "$(dirname "$0")/.."
JULIA=${JULIA:-julia}
JULIA_VERSION=$($JULIA --version 2>/dev/null | awk 'NR==1 {print $3}' || true)
if [[ $JULIA_VERSION != 1.6.* ]]; then
  echo "This project needs Julia 1.6 (found '${JULIA_VERSION:-no julia}' from '$JULIA')."
  echo "With juliaup: juliaup add 1.6, then run with JULIA=\"julia +1.6\" $0 ..."
  exit 1
fi

OUTDIR="results/$(date +%Y-%m-%d)"
REPS=100000
SENS_REPS=10000
SEASONS=all
GAMMA=auto
MODE=-2
JOBS=1
THREADS=
EXPERIMENTS=validate,simulate,parse,noisy,sensitivity
PLOT=--plot

while getopts "o:n:S:s:g:m:j:t:e:Nh" opt; do
  case $opt in
    o) OUTDIR=$OPTARG ;;
    n) REPS=$OPTARG ;;
    S) SENS_REPS=$OPTARG ;;
    s) SEASONS=$OPTARG ;;
    g) GAMMA=$OPTARG ;;
    m) MODE=$OPTARG ;;
    j) JOBS=$OPTARG ;;
    t) THREADS=$OPTARG ;;
    e) EXPERIMENTS=$OPTARG ;;
    N) PLOT= ;;
    h) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 0 ;;
    *) sed -n '2,/^set -euo/p' "$0" | sed '$d'; exit 1 ;;
  esac
done

[[ -n $THREADS ]] || THREADS=$(( JOBS > 1 ? JOBS : 4 ))
if (( JOBS > 1 )) && [[ -z ${TANKING_GUROBI_THREADS:-} ]]; then
  export TANKING_GUROBI_THREADS=1
fi
has() { [[ ",$EXPERIMENTS," == *",$1,"* ]]; }
RUN="$JULIA --project=. scripts/run_experiments.jl --results-dir=$OUTDIR --replications=$REPS --seasons=$SEASONS"
LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"
NUM_STEPS=31 # 0, 1, ..., 30 selfish teams

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }
# show_error LOGFILE: show the (first) error message in a log, or its last lines if there is none
show_error() {
  if grep -q "ERROR" "$1"; then
    grep -m1 -A4 "ERROR" "$1" | sed 's/^/    /'
  else
    tail -n 10 "$1" | sed 's/^/    /'
  fi
}
# run_parallel COUNT LOGPREFIX COMMAND: run COMMAND for {} = 1..COUNT in $JOBS parallel jobs, with output
# to LOGPREFIX{}.log; on failure, show the error from a failed job's log and stop
run_parallel() {
  local count=$1 logprefix=$2 cmd=$3
  if ! seq 1 "$count" | xargs -P "$JOBS" -I{} sh -c \
      "$cmd > ${logprefix}{}.log 2>&1 || { echo 'job {} failed; see ${logprefix}{}.log'; exit 255; }"; then
    failed=$(grep -l "ERROR" "${logprefix}"*.log | head -1)
    log "FAILED; see ${failed:-the logs ${logprefix}*.log}:"
    [[ -n $failed ]] && show_error "$failed"
    exit 1
  fi
}
# run_step LOGFILE COMMAND...: run COMMAND with output to LOGFILE; on failure, show the end of the log and stop
run_step() {
  local logfile=$1; shift
  if ! "$@" > "$logfile" 2>&1; then
    log "FAILED; see $logfile:"
    show_error "$logfile"
    exit 1
  fi
}

{
  echo "date: $(date)"
  echo "git commit: $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "julia: $JULIA_VERSION ($JULIA)"
  echo "replications: $REPS (sensitivity: $SENS_REPS)"
  echo "seasons: $SEASONS"
  echo "gamma: $GAMMA"
  echo "math_elim_mode: $MODE"
  echo "jobs: $JOBS, threads: $THREADS, Gurobi threads per MIP: ${TANKING_GUROBI_THREADS:-default}"
  echo "CPU affinity: $(taskset -pc $$ 2>/dev/null | sed 's/.*: //' || echo unknown)"
  echo "experiments: $EXPERIMENTS"
  echo "plots: ${PLOT:-no}"
} > "$OUTDIR/settings.txt"
log "Settings:"; sed 's/^/    /' "$OUTDIR/settings.txt"

## 1. Model validation (choice of gamma)
if has validate; then
  GAMMA_OPT=""
  [[ $GAMMA == auto ]] && GAMMA_OPT="--gamma=auto"
  if (( JOBS > 1 )); then
    NUM_MODELS=$($JULIA --project=. -e 'using Tanking; print(Tanking.num_validation_models())' 2>/dev/null | tail -1)
    [[ $NUM_MODELS =~ ^[0-9]+$ ]] || { log "FAILED to get the number of validation models (got: $NUM_MODELS)"; exit 1; }
    log "validate: $NUM_MODELS models in $JOBS parallel jobs -> $LOGDIR/validate_model*.log"
    run_parallel "$NUM_MODELS" "$LOGDIR/validate_model" "$RUN --model={} validate"
    log "validate: combining models -> $LOGDIR/validate.log"
    run_step "$LOGDIR/validate.log" $RUN $GAMMA_OPT --aggregate $PLOT validate
  else
    log "validate: model_validation -> $LOGDIR/validate.log"
    run_step "$LOGDIR/validate.log" $RUN $GAMMA_OPT $PLOT validate
  fi
  log "validate: minimax gamma = $(cat "$OUTDIR/gamma.txt")"
  if grep -q "Warning:" "$LOGDIR/validate.log"; then
    grep -o "Warning: .*" "$LOGDIR/validate.log" | sed 's/^/    /'
  fi
fi
if [[ $GAMMA == auto ]]; then
  [[ -f $OUTDIR/gamma.txt ]] || { echo "gamma=auto needs the validate experiment (or $OUTDIR/gamma.txt)"; exit 1; }
  GAMMA=$(cat "$OUTDIR/gamma.txt")
  echo "gamma (auto): $GAMMA" >> "$OUTDIR/settings.txt"
fi

## 2. Simulation (STRICT mode), and 6. Bradley-Terry simulation, possibly split over steps
simulate_experiment() { # $1 = simulate or bt
  local exp=$1
  if (( JOBS > 1 )); then
    log "$exp: $NUM_STEPS steps in $JOBS parallel jobs -> $LOGDIR/${exp}_step*.log"
    run_parallel "$NUM_STEPS" "$LOGDIR/${exp}_step" "$RUN --gamma=$GAMMA --math-elim-mode=$MODE --steps={} $exp"
    log "$exp: aggregating steps -> $LOGDIR/${exp}_aggregate.log"
    run_step "$LOGDIR/${exp}_aggregate.log" $RUN --gamma=$GAMMA --math-elim-mode=$MODE --aggregate $PLOT $exp
  else
    log "$exp: all steps in one process -> $LOGDIR/$exp.log"
    run_step "$LOGDIR/$exp.log" $RUN --gamma=$GAMMA --math-elim-mode=$MODE $PLOT $exp
  fi
}
if has simulate; then simulate_experiment simulate; fi

## 3. NBA data (needs avg_eff_eliminated_strict.csv from the simulation)
if has parse; then
  log "parse: main_parse -> $LOGDIR/parse.log"
  run_step "$LOGDIR/parse.log" $RUN $PLOT parse
fi

## 4. Noisy rankings
if has noisy; then
  log "noisy: rankings_are_noisy -> $LOGDIR/noisy.log"
  run_step "$LOGDIR/noisy.log" $RUN $PLOT noisy
fi

## 5. Sensitivity to tanking after the breakpoint (effective elimination; no Gurobi needed)
if has sensitivity; then
  log "sensitivity: $SENS_REPS replications, $THREADS threads -> $LOGDIR/sensitivity.log"
  run_step "$LOGDIR/sensitivity.log" $JULIA --project=. -t "$THREADS" scripts/run_sensitivity.jl "$SENS_REPS" "$OUTDIR/sensitivity" "1:$NUM_STEPS" 628 "$GAMMA"
  grep "^Check" "$LOGDIR/sensitivity.log" | sed 's/^/    /'
fi

## 6. Bradley-Terry model estimated from the NBA seasons
if has bt; then simulate_experiment bt; fi

log "Done. Results in $OUTDIR"
