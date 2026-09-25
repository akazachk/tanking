#!/usr/bin/env julia
#
# Check the output of scripts/run_all.sh (used by test/test_all.sh):
#   julia --project=. test/check_results.jl <results_dir>
#
# Checks that every expected file exists, has the right size, and contains finite values;
# that min <= avg <= max and stddev >= 0 for every statistic; that some basic properties hold
# (e.g., no games are tanked without selfish teams, games tanked do not decrease with the breakpoint);
# that Gurobi was actually used when math_elim_mode requires it (MIPs were solved);
# that the requested NBA seasons were used; and that the sensitivity checks are 0.
# Exits with status 1 if any check fails.

using DelimitedFiles
using Printf

const NUM_STEPS = 31        # 0, ..., 30 selfish teams
const NUM_BREAKPOINTS = 6
const NUM_NBA_ODDS = 3
const NUM_GAMES = 3 * 30 * 29 ÷ 2
const PREFIX = ["avg_", "stddev_", "min_", "max_"]
const TOL = 1e-6

num_fail = 0
num_pass = 0
function check(ok::Bool, msg)
  global num_fail, num_pass
  if ok
    num_pass += 1
  else
    num_fail += 1
    println("FAIL: ", msg)
  end
  return ok
end

function read_matrix(file)
  M = readdlm(file, ',')
  return M isa AbstractVector ? reshape(M, :, 1) : M
end

function main(dir)
  settings = Dict{String,String}()
  settings_file = joinpath(dir, "settings.txt")
  if !check(isfile(settings_file), "$settings_file not found (was scripts/run_all.sh run?)")
    return
  end
  for line in eachline(settings_file)
    k, v = strip.(split(line, ":", limit=2))
    settings[k] = v
  end
  mode = parse(Int, settings["math_elim_mode"])
  experiments = split(settings["experiments"], ",")
  plots = settings["plots"] != "no"
  seasons = settings["seasons"]
  println("Checking $dir (math_elim_mode = $mode, seasons = $seasons, experiments = $(join(experiments, ",")), plots = $plots)")

  ## Simulation output
  if "simulate" in experiments
    sizes = Dict("kend" => (NUM_STEPS, NUM_BREAKPOINTS), "kend_nba" => (NUM_STEPS, NUM_NBA_ODDS),
        "kend_lenten" => (NUM_STEPS, 1), "games_tanked" => (NUM_STEPS, NUM_BREAKPOINTS),
        "already_tank" => (NUM_STEPS, NUM_BREAKPOINTS), "math_eliminated" => (NUM_STEPS, NUM_GAMES),
        "eff_eliminated" => (NUM_STEPS, NUM_GAMES), "num_mips" => (NUM_STEPS, 1), "num_unelim" => (NUM_STEPS, 1),
        "avg_rank_strat" => (NUM_STEPS, 1), "avg_rank_moral" => (NUM_STEPS, 1),
        "avg_elim_rank_strat" => (NUM_STEPS, 1), "avg_elim_rank_moral" => (NUM_STEPS, 1),
        "avg_diff_rank_strat" => (NUM_STEPS, 1), "avg_diff_rank_moral" => (NUM_STEPS, 1),
        "num_missing_case" => (NUM_STEPS, 1))
    data = Dict{String,Vector{Matrix{Float64}}}()
    for (name, sz) in sizes
      mats = Matrix{Float64}[]
      for p in PREFIX
        file = joinpath(dir, string(p, name, "_strict.csv"))
        if !check(isfile(file), "missing $file")
          continue
        end
        M = Float64.(read_matrix(file))
        check(size(M) == sz, "$file has size $(size(M)), expected $sz")
        check(all(isfinite, M), "$file has non-finite values")
        push!(mats, M)
      end
      length(mats) == 4 || continue
      data[name] = mats
      avg, sd, mn, mx = mats
      # Statistics that are only averaged over some replications may keep their initial min/max when never updated
      partial = name in ("avg_elim_rank_strat", "avg_elim_rank_moral", "avg_diff_rank_strat", "avg_diff_rank_moral")
      if !partial
        check(all(mn .<= avg .+ TOL) && all(avg .<= mx .+ TOL), "$name: min <= avg <= max violated")
      end
      check(all(sd .>= -TOL), "$name: negative standard deviation")
    end

    if haskey(data, "games_tanked")
      gt = data["games_tanked"][1]
      check(all(gt[1, :] .== 0), "games tanked with 0 selfish teams should be 0, got $(gt[1, :])")
      check(all(diff(gt, dims=2) .>= -TOL), "games tanked should not decrease with the breakpoint")
      check(maximum(gt) > 0, "no games tanked in any step")
    end
    if haskey(data, "kend")
      k = data["kend"][1]
      check(all(0 .< k .< 91), "average Kendall tau outside (0, 91)")
    end
    if haskey(data, "eff_eliminated")
      e = data["eff_eliminated"][1][:, end]
      check(all(0 .<= e .<= 30), "number effectively eliminated at end of season outside [0, 30]")
      check(maximum(e) > 0, "no teams effectively eliminated")
    end
    if haskey(data, "kend_lenten")
      l = data["kend_lenten"][1]
      check(all(0 .< l .< 91), "average Lenten Kendall tau outside (0, 91)")
    end
    gold_file = joinpath(dir, "kend_gold_strict.csv")
    if check(isfile(gold_file), "missing $gold_file")
      g = read_matrix(gold_file)
      check(length(g) == 4 && g[1] > 0, "Gold Kendall tau should be positive (is it all zeros after aggregation?): $g")
    end
    if abs(mode) >= 2 && haskey(data, "num_mips")
      check(maximum(data["num_mips"][1]) > 0, "math_elim_mode = $mode but no MIPs were solved (was Gurobi used?)")
    end
    if mode != 0 && haskey(data, "math_eliminated")
      check(maximum(data["math_eliminated"][1][:, end]) > 0, "math_elim_mode = $mode but no teams were mathematically eliminated")
      check(maximum(data["num_missing_case"][1]) >= 0, "num_missing_case negative")
    end
    if mode == 0 && haskey(data, "num_mips")
      check(maximum(data["num_mips"][1]) == 0, "math_elim_mode = 0 but MIPs were solved")
    end
    if haskey(data, "kend")
      avg = data["kend"][1]; avg_n = data["num_missing_case"][1]; mips = data["num_mips"][1]
      @printf("  simulate: avg Kendall tau (0 selfish, end of season) = %.2f; avg num_missing_case = %.2f; avg MIPs per replication = %.1f\n",
          avg[1, end], sum(avg_n) / length(avg_n), sum(mips) / length(mips))
    end
    if plots
      for f in ["avg_kend", "avg_games_tanked", "avg_eliminated", "avg_rank", "avg_elim_rank", "avg_already_tank"]
        check(isfile(joinpath(dir, "pdf", f * "_strict.pdf")), "missing plot pdf/$(f)_strict.pdf")
      end
    end
  end

  ## Model validation
  if "validate" in experiments
    for f in ["model_validation_strict.csv", "gamma.txt"]
      check(isfile(joinpath(dir, f)), "missing $(joinpath(dir, f))")
    end
    if haskey(settings, "gamma (auto)") && isfile(joinpath(dir, "gamma.txt"))
      g = strip(read(joinpath(dir, "gamma.txt"), String))
      check(settings["gamma (auto)"] == g, "gamma used ($(settings["gamma (auto)"])) differs from gamma.txt ($g)")
      println("  validate: gamma chosen by model validation = $g")
    end
    log = joinpath(dir, "logs", "validate.log")
    if check(isfile(log), "missing $log")
      text = read(log, String)
      expected = seasons == "all" ? "2025-26" : "2018-19"
      check(occursin(expected, text), "validate.log does not mention season $expected")
      if seasons != "all"
        check(!occursin("2021-22", text), "validate.log mentions 2021-22 although seasons = $seasons")
      end
    end
    if plots
      for f in ["model_loss", "win_pct_0tank"]
        check(isfile(joinpath(dir, "pdf", f * "_strict.pdf")), "missing plot pdf/$(f)_strict.pdf")
      end
    end
  end

  ## NBA data
  if "parse" in experiments
    log = joinpath(dir, "logs", "parse.log")
    check(isfile(log) && !occursin("ERROR", read(log, String)), "parse failed; see $log")
    if plots
      for f in ["nba_num_games_tanked", "nba_num_teams_eliminated"]
        check(isfile(joinpath(dir, "pdf", f * "_strict.pdf")), "missing plot pdf/$(f)_strict.pdf")
      end
    end
  end

  ## Noisy rankings
  if "noisy" in experiments
    file = joinpath(dir, "noisy_ranking_strict.csv")
    if check(isfile(file), "missing $file")
      M = read_matrix(file)
      check(size(M) == (51, 8) && all(isfinite, M), "$file has size $(size(M)) or non-finite values")
    end
    plots && check(isfile(joinpath(dir, "pdf", "noisy_ranking_strict.pdf")), "missing plot pdf/noisy_ranking_strict.pdf")
  end

  ## Sensitivity
  if "sensitivity" in experiments
    file = joinpath(dir, "sensitivity", "checks.txt")
    if check(isfile(file), "missing $file")
      for line in eachline(file)
        if startswith(line, "max")
          val = parse(Float64, strip(split(line, "=")[end]))
          check(val == 0, "sensitivity check not 0: $line")
        end
      end
    end
    for f in ["kend_keep.csv", "kend_stop.csv", "kend_diff.csv", "se_diff.csv"]
      check(isfile(joinpath(dir, "sensitivity", f)), "missing sensitivity/$f")
    end
  end

  ## Bradley-Terry
  if "bt" in experiments
    check(isfile(joinpath(dir, "avg_kend_BT_est.csv")), "missing avg_kend_BT_est.csv")
  end

  println(num_fail == 0 ? "PASSED" : "FAILED", ": $num_pass checks passed, $num_fail failed")
end

if abspath(PROGRAM_FILE) == @__FILE__
  if length(ARGS) != 1
    println("Usage: julia --project=. test/check_results.jl <results_dir>")
    exit(2)
  end
  main(ARGS[1])
  exit(num_fail == 0 ? 0 : 1)
end
