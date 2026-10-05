#!/usr/bin/env julia
#
# Sensitivity of the Kendall tau results to tanking after the breakpoint.
#
# In `simulate`, a selfish team tanks in every game after it is (effectively) eliminated,
# including games after the breakpoint delta, so that one season serves all breakpoints.
# Here we compare, on identical seasons up to delta (common random numbers):
#   keep: the default behavior (one run over all breakpoints)
#   stop: no team tanks after delta (one run per breakpoint)
# using STRICT mode and effective elimination (math_elim_mode = 0; tanking decisions are
# the same as with the default math_elim_mode = -2, and Gurobi is not needed).
#
# Usage (from the project root):
#   julia --project=. [-t <threads>] scripts/run_sensitivity.jl <num_replications> <results_dir> [steps] [seed] [gamma]
# e.g.
#   julia --project=. -t 4 scripts/run_sensitivity.jl 10000 results/sens_test "[1,9,16,24,31]"
# Step s corresponds to s-1 selfish teams (default: all steps 1:31). Default seed: 628.
# gamma: if not given, it is read from gamma.txt (written by model validation) in <results_dir> or its parent
# directory (e.g., results/<date>/gamma.txt for results_dir = results/<date>/sensitivity).
#
# Output (rows: steps, columns: breakpoints; first column is the number of selfish teams):
#   kend_keep.csv, kend_stop.csv, kend_diff.csv (= stop - keep): average Kendall tau distance
#   se_keep.csv, se_stop.csv, se_diff.csv: standard errors of the above (se_diff is from the paired differences)
#   games_tanked_keep.csv, games_tanked_stop.csv: average number of games tanked up to each breakpoint
#   breakpoints.csv: the breakpoints (fraction of season) and corresponding games
# Checks (printed; both must be exactly 0):
#   max |games tanked (keep) - games tanked (stop)| over all replications, steps, and breakpoints
#   max |Kendall tau (keep) - Kendall tau (stop)| at delta = end of season

using Tanking
using DelimitedFiles
using Printf

function str2range(input::AbstractString)
  if findfirst(':', input) != nothing
    y = split(input, ':')
    return UnitRange(parse.(Int,y[1]),parse.(Int,y[2]))
  else
    return parse.(Int,input)
  end
end # str2range

function str2arr(input::AbstractString)
  y = split(input, x -> (x == '[' || x == ',' || x == ']' || isspace(x)) ? true : false)
  y = [y[i] for i in 1:length(y) if y[i] != ""]
  return vcat([collect(str2range(y[i])) for i in 1:length(y)]...)
end # str2arr

"""
Run one replication (seeded) of `simulate`, returning the Kendall tau distance and number of games tanked
for each breakpoint in `bp_list`
"""
function one_replication(step, seed, bp_list, stop, gamma; num_rounds=3)
  out = Tanking.simulate(Tanking.num_teams, Tanking.num_playoff_teams, num_rounds, 1, Tanking.num_teams, gamma,
      bp_list, Tanking.nba_odds_list, Tanking.nba_num_lottery, Tanking.true_strength, Tanking.STRICT,
      0, [step], nothing, false; stop_tanking_after_breakpoint=stop, seed_per_replication=seed, verbose=false)
  kend, games_tanked = out[1], out[5]
  avg_stat = 1
  return kend[step, :, avg_stat], games_tanked[step, :, avg_stat]
end # one_replication

function main(args)
  num_replications = parse(Int, args[1])
  results_dir = args[2]
  steps = length(args) >= 3 ? str2arr(args[3]) : collect(1:Tanking.num_teams+1)
  seed = length(args) >= 4 ? parse(Int, args[4]) : 628
  gamma = nothing
  if length(args) >= 5
    gamma = parse(Float64, args[5])
  else
    for dir in (results_dir, dirname(normpath(results_dir)))
      file = joinpath(dir, "gamma.txt")
      if isfile(file)
        gamma = parse(Float64, strip(read(file, String)))
        println("Using gamma = $gamma from $file")
        break
      end
    end
    isnothing(gamma) && error("No gamma given and no gamma.txt found in $results_dir or its parent; pass gamma as the 5th argument")
  end
  mkpath(results_dir)
  Tanking.set_mode(Tanking.STRICT)

  bp = Tanking.breakpoint_list
  num_bp = length(bp)
  num_steps = length(steps)
  println("Sensitivity run: $num_replications replications, steps $steps, breakpoints $bp, seed $seed, gamma $gamma, $(Threads.nthreads()) threads")

  # Per-replication values [rep, step, breakpoint]
  kend_keep = zeros(num_replications, num_steps, num_bp)
  kend_stop = zeros(num_replications, num_steps, num_bp)
  tanked_keep = zeros(num_replications, num_steps, num_bp)
  tanked_stop = zeros(num_replications, num_steps, num_bp)

  for (s_ind, step) in enumerate(steps)
    t = @elapsed Threads.@threads for rep = 1:num_replications
      # Replication rep of this step uses seed + rep in every run (simulate adds 1 to the seed for its only replication)
      curr_seed = seed + 1_000_000 * step + rep - 1
      k, g = one_replication(step, curr_seed, bp, false, gamma)
      kend_keep[rep, s_ind, :] = k
      tanked_keep[rep, s_ind, :] = g
      for r = 1:num_bp
        k, g = one_replication(step, curr_seed, [bp[r]], true, gamma)
        kend_stop[rep, s_ind, r] = k[1]
        tanked_stop[rep, s_ind, r] = g[1]
      end
    end
    @printf("Step %d (%d selfish teams) done in %.1f s\n", step, step-1, t)
  end

  ## Checks
  check_tanked = maximum(abs.(tanked_keep - tanked_stop))
  end_bp = findfirst(isequal(1), bp)
  check_end = isnothing(end_bp) ? NaN : maximum(abs.(kend_keep[:, :, end_bp] - kend_stop[:, :, end_bp]))
  println("Check: max |games tanked (keep) - games tanked (stop)| = ", check_tanked, " (must be 0)")
  println("Check: max |Kendall tau (keep) - Kendall tau (stop)| at delta = T: ", check_end, " (must be 0)")

  ## Summaries
  mean_over_reps(x) = dropdims(sum(x, dims=1), dims=1) / num_replications
  function se_over_reps(x)
    m = mean_over_reps(x)
    v = dropdims(sum((x .- reshape(m, 1, size(m)...)).^2, dims=1), dims=1) / max(num_replications - 1, 1)
    return sqrt.(v / num_replications)
  end
  diff = kend_stop - kend_keep

  header = hcat("num_selfish", reshape([string(b) for b in bp], 1, num_bp))
  function write_table(name, M)
    writedlm(joinpath(results_dir, name), vcat(header, hcat(steps .- 1, M)), ',')
  end
  write_table("kend_keep.csv", mean_over_reps(kend_keep))
  write_table("kend_stop.csv", mean_over_reps(kend_stop))
  write_table("kend_diff.csv", mean_over_reps(diff))
  write_table("se_keep.csv", se_over_reps(kend_keep))
  write_table("se_stop.csv", se_over_reps(kend_stop))
  write_table("se_diff.csv", se_over_reps(diff))
  write_table("games_tanked_keep.csv", mean_over_reps(tanked_keep))
  write_table("games_tanked_stop.csv", mean_over_reps(tanked_stop))
  num_games_total = 3 * Tanking.num_teams * (Tanking.num_teams - 1) ÷ 2
  writedlm(joinpath(results_dir, "breakpoints.csv"),
      vcat(["breakpoint" "game"], hcat(string.(bp), [round(b * num_games_total) for b in bp])), ',')
  open(joinpath(results_dir, "checks.txt"), "w") do io
    println(io, "replications = $num_replications, seed = $seed, gamma = $gamma, steps = $steps")
    println(io, "max |games tanked (keep) - games tanked (stop)| = $check_tanked")
    println(io, "max |Kendall tau (keep) - Kendall tau (stop)| at delta = T = $check_end")
  end

  println("\nAverage change in Kendall tau (stop - keep), with standard errors:")
  md = mean_over_reps(diff); sd = se_over_reps(diff)
  @printf("%12s", "selfish")
  for b in bp
    @printf("%18s", string(b))
  end
  println()
  for (s_ind, step) in enumerate(steps)
    @printf("%12d", step - 1)
    for r = 1:num_bp
      @printf("%18s", @sprintf("%.3f (%.3f)", md[s_ind, r], sd[s_ind, r]))
    end
    println()
  end
  println("\nResults written to $results_dir")
end # main

if abspath(PROGRAM_FILE) == @__FILE__
  if length(ARGS) < 2
    println("Usage: julia --project=. [-t threads] scripts/run_sensitivity.jl <num_replications> <results_dir> [steps] [seed] [gamma]")
  else
    main(ARGS)
  end
end
