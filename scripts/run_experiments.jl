#!/usr/bin/env julia
#
# Rerun the experiments of the paper, using the NBA seasons listed in
# `Tanking.nba_seasons` (2004-05 to 2025-26, except 2011-12, 2019-20, 2020-21).
#
# Usage (from the project root):
#   julia --project=. scripts/run_experiments.jl [options] [experiments...]
#
# Experiments (run in this order; default: all of them)
#   validate  model_validation: compare simulated win pct by rank (STRICT mode with several
#             values of gamma, and Bradley-Terry MLE) to NBA data in data/winpct.csv
#   simulate  main_simulate: tanking simulation (STRICT mode)
#   parse     main_parse: games that could be tanked in real NBA seasons
#             (needs avg_eff_eliminated_strict.csv from `simulate` in the same results directory)
#   noisy     rankings_are_noisy: noisiness of the ranking as a function of gamma and number of rounds
#   bt        main_simulate with the Bradley-Terry model estimated from NBA data (BT_ESTIMATED mode)
#             (not part of the default set)
#
# Options
#   --results-dir=DIR       where results are written (default: results/<today's date>)
#   --replications=N        number of replications (default: 100000)
#   --gamma=G               probability the better team wins in STRICT mode (default: 0.71425);
#                           use --gamma=auto to use the value whose largest model_validation loss (over 0, 15, 30
#                           selfish teams) is smallest, i.e., the minimax rule used for 0.71425
#                           (requires `validate` to be run first, in the same call)
#   --math-elim-mode=M      see README (default: -2)
#   --steps=S               only simulate step S (e.g., 5); to split the simulation across jobs, run
#                           once for each step (1 to 31), then once with --aggregate (only single-step
#                           runs can be aggregated)
#   --aggregate             combine results of runs that used --steps (main_simulate with do_simulation=-1)
#   --plot                  also create plots (needs PyPlot and LaTeX)
#   --seasons=S             NBA seasons used by validate, parse, and bt: "all" (default; 2004-05 to 2025-26,
#                           except 2011-12, 2019-20, 2020-21) or "2004-2019" (the 14 seasons of the original
#                           paper, 2004-05 to 2018-19); the simulations themselves do not use NBA data
#
# Examples
#   julia --project=. scripts/run_experiments.jl --replications=100 --results-dir=results/tmp
#   julia --project=. scripts/run_experiments.jl --gamma=auto validate simulate parse

using Tanking
using Dates

"""
    str2range

Parse string (input should be a single number or range as a string) into UnitRange or into just an Int if it is contains no colon
"""
function str2range(input::AbstractString)
  if findfirst(':', input) != nothing
    y = split(input, ':')
    return UnitRange(parse.(Int,y[1]),parse.(Int,y[2]))
  else
    return parse.(Int,input)
  end
end # str2range

"""
    str2arr

Convert string into array
"""
function str2arr(input::AbstractString)
  y = split(input, x -> (x == '[' || x == ',' || x == ']' || isspace(x)) ? true : false)
  y = [y[i] for i in 1:length(y) if y[i] != ""]
  return [str2range(y[i]) for i in 1:length(y)]
end # str2arr

const ALL_EXPERIMENTS = ["validate", "simulate", "parse", "noisy", "bt"]
const DEFAULT_EXPERIMENTS = ["validate", "simulate", "parse", "noisy"]

function main(args)
  results_dir = joinpath("results", string(Dates.today()))
  num_replications = 100000
  gamma_arg = "0.71425"
  math_elim_mode = -2
  selected_steps = nothing
  aggregate = false
  do_plotting = false
  experiments = String[]

  for arg in args
    if startswith(arg, "--results-dir=")
      results_dir = split(arg, "=", limit=2)[2]
    elseif startswith(arg, "--replications=")
      num_replications = parse(Int, split(arg, "=", limit=2)[2])
    elseif startswith(arg, "--gamma=")
      gamma_arg = split(arg, "=", limit=2)[2]
    elseif startswith(arg, "--math-elim-mode=")
      math_elim_mode = parse(Int, split(arg, "=", limit=2)[2])
    elseif startswith(arg, "--steps=")
      selected_steps = str2arr(split(arg, "=", limit=2)[2])
    elseif arg == "--aggregate"
      aggregate = true
    elseif startswith(arg, "--seasons=")
      val = split(arg, "=", limit=2)[2]
      if val == "all"
        Tanking.set_seasons!(Tanking.nba_seasons)
      elseif val == "2004-2019"
        Tanking.set_seasons!(Tanking.nba_seasons_2004_2019)
      else
        error("Unknown value for --seasons: $val (use all or 2004-2019)")
      end
    elseif arg == "--plot"
      do_plotting = true
    elseif arg in ALL_EXPERIMENTS
      push!(experiments, arg)
    else
      error("Unknown argument: $arg")
    end
  end
  if isempty(experiments)
    experiments = DEFAULT_EXPERIMENTS
  end
  experiments = [e for e in ALL_EXPERIMENTS if e in experiments] # run in canonical order
  if gamma_arg == "auto" && !("validate" in experiments)
    error("--gamma=auto requires the validate experiment")
  end
  gamma = gamma_arg == "auto" ? nothing : parse(Float64, gamma_arg)

  mkpath(results_dir)
  println("## Running experiments ", experiments, " with results in ", results_dir)
  println("## NBA seasons: ", Tanking.selected_nba_seasons)

  if "validate" in experiments
    println("\n## model_validation ##")
    # best_gamma: chosen by the minimax rule (smallest largest loss over 0, 15, 30 selfish teams),
    # as used to choose 0.71425 in the paper; it is also the value drawn in the win-pct plots
    @time loss_list, gamma_list, best_gamma = Tanking.model_validation(do_simulation=true,
        num_replications=num_replications, results_dir=results_dir,
        do_plotting=do_plotting, selected_steps=nothing)
    println("Value of gamma with smallest maximum loss over 0, 15, 30 selfish teams: ", best_gamma)
    open(joinpath(results_dir, "gamma.txt"), "w") do io
      println(io, best_gamma)
    end
    if isnothing(gamma)
      gamma = best_gamma
    end
  end
  if !isnothing(gamma)
    println("## Using gamma = ", gamma)
  end

  if "simulate" in experiments
    println("\n## main_simulate (STRICT) ##")
    @time Tanking.main_simulate(do_simulation=(aggregate ? -1 : 1),
        num_replications=num_replications, do_plotting=do_plotting,
        mode=Tanking.STRICT, results_dir=results_dir, gamma=gamma,
        math_elim_mode=math_elim_mode, selected_steps=(aggregate ? nothing : selected_steps))
  end

  if "parse" in experiments
    println("\n## main_parse ##")
    @time Tanking.main_parse(do_plotting=do_plotting, mode=Tanking.STRICT, results_dir=results_dir)
  end

  if "noisy" in experiments
    println("\n## rankings_are_noisy ##")
    @time Tanking.rankings_are_noisy(do_simulation=true, num_replications=num_replications,
        do_plotting=do_plotting, mode=Tanking.STRICT, results_dir=results_dir)
  end

  if "bt" in experiments
    println("\n## main_simulate (BT_ESTIMATED) ##")
    @time Tanking.main_simulate(do_simulation=(aggregate ? -1 : 1),
        num_replications=num_replications, do_plotting=do_plotting,
        mode=Tanking.BT_ESTIMATED, results_dir=results_dir,
        math_elim_mode=math_elim_mode, selected_steps=(aggregate ? nothing : selected_steps))
  end
end # main

if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS)
end
