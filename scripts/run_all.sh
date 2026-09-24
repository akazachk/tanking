#!/usr/bin/env bash
#
# Recreate all experiments and plots, from the main project directory:
#
#   scripts/run_all.sh [options]
#
# Options
#   -o DIR     results directory (default: results/run_<today>)
#   -n N       replications for model validation, the simulation, and noisy rankings (default: 100000)
#   -S N       replications for the sensitivity run (default: 10000)
#   -s SEASONS NBA seasons for model validation, NBA parsing, and Bradley-Terry: all (default) or 2004-2019
#   -g GAMMA   gamma for the simulation (default: 0.71425), or "auto" for the minimax choice from model validation
#   -m MODE    math_elim_mode for the simulation (default: -2, which needs Gurobi)
#   -j JOBS    number of parallel jobs for the simulation (default: 1); with JOBS > 1, each of the
#              31 steps is simulated in its own process and the results are then aggregated
#   -t THREADS threads for the sensitivity run (default: 4)
#   -e LIST    comma-separated experiments to run, from
#              validate,simulate,parse,noisy,sensitivity,bt (default: validate,simulate,parse,noisy,sensitivity)
#   -N         do not create plots (plots need PyPlot and LaTeX)
#
# Environment variables
#   JULIA      julia command (default: julia); e.g., JULIA="julia --sysimage=build/JuliaTanking.so"
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

OUTDIR="results/run_$(date +%Y-%m-%d)"
REPS=100000
SENS_REPS=10000
SEASONS=all
GAMMA=0.71425
MODE=-2
JOBS=1
THREADS=4
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
    h) sed -n '2,32p' "$0"; exit 0 ;;
    *) sed -n '2,32p' "$0"; exit 1 ;;
  esac
done

has() { [[ ",$EXPERIMENTS," == *",$1,"* ]]; }
RUN="$JULIA --project=. scripts/run_experiments.jl --results-dir=$OUTDIR --replications=$REPS --seasons=$SEASONS"
LOGDIR="$OUTDIR/logs"
mkdir -p "$LOGDIR"
NUM_STEPS=31 # 0, 1, ..., 30 selfish teams

log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"; }

{
  echo "date: $(date)"
  echo "git commit: $(git rev-parse HEAD 2>/dev/null || echo unknown)"
  echo "replications: $REPS (sensitivity: $SENS_REPS)"
  echo "seasons: $SEASONS"
  echo "gamma: $GAMMA"
  echo "math_elim_mode: $MODE"
  echo "jobs: $JOBS, threads: $THREADS"
  echo "experiments: $EXPERIMENTS"
  echo "plots: ${PLOT:-no}"
} > "$OUTDIR/settings.txt"
log "Settings:"; sed 's/^/    /' "$OUTDIR/settings.txt"

## 1. Model validation (choice of gamma)
if has validate; then
  log "validate: model_validation -> $LOGDIR/validate.log"
  GAMMA_OPT=""
  [[ $GAMMA == auto ]] && GAMMA_OPT="--gamma=auto"
  $RUN $GAMMA_OPT $PLOT validate > "$LOGDIR/validate.log" 2>&1
  log "validate: minimax gamma = $(cat "$OUTDIR/gamma.txt")"
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
    seq 1 $NUM_STEPS | xargs -P "$JOBS" -I{} sh -c \
      "$RUN --gamma=$GAMMA --math-elim-mode=$MODE --steps={} $exp > $LOGDIR/${exp}_step{}.log 2>&1 || { echo 'step {} failed; see $LOGDIR/${exp}_step{}.log'; exit 255; }"
    log "$exp: aggregating steps -> $LOGDIR/${exp}_aggregate.log"
    $RUN --gamma=$GAMMA --math-elim-mode=$MODE --aggregate $PLOT $exp > "$LOGDIR/${exp}_aggregate.log" 2>&1
  else
    log "$exp: all steps in one process -> $LOGDIR/$exp.log"
    $RUN --gamma=$GAMMA --math-elim-mode=$MODE $PLOT $exp > "$LOGDIR/$exp.log" 2>&1
  fi
}
if has simulate; then simulate_experiment simulate; fi

## 3. NBA data (needs avg_eff_eliminated_strict.csv from the simulation)
if has parse; then
  log "parse: main_parse -> $LOGDIR/parse.log"
  $RUN $PLOT parse > "$LOGDIR/parse.log" 2>&1
fi

## 4. Noisy rankings
if has noisy; then
  log "noisy: rankings_are_noisy -> $LOGDIR/noisy.log"
  $RUN $PLOT noisy > "$LOGDIR/noisy.log" 2>&1
fi

## 5. Sensitivity to tanking after the breakpoint (effective elimination; no Gurobi needed)
if has sensitivity; then
  log "sensitivity: $SENS_REPS replications, $THREADS threads -> $LOGDIR/sensitivity.log"
  $JULIA --project=. -t "$THREADS" scripts/run_sensitivity.jl "$SENS_REPS" "$OUTDIR/sensitivity" "1:$NUM_STEPS" 628 "$GAMMA" > "$LOGDIR/sensitivity.log" 2>&1
  grep "^Check" "$LOGDIR/sensitivity.log" | sed 's/^/    /'
fi

## 6. Bradley-Terry model estimated from the NBA seasons
if has bt; then simulate_experiment bt; fi

log "Done. Results in $OUTDIR"
