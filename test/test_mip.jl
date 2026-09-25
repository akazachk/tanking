#!/usr/bin/env julia
#
# Test the mathematical-elimination MIPs without a Gurobi license, using the open-source solver HiGHS:
#   julia --project=test test/test_mip.jl [num_replications]
# (from the main project directory; the test environment in test/Project.toml adds HiGHS and uses the local Tanking)
#
# Simulates seasons with math_elim_mode = 2 (binary team-wise MIP) while
#   * DEBUG = true: every elimination check is also done with the general-integer formulation (mode 3),
#     and the two exact formulations must agree (asserted in simulate);
#   * CHECK_BEST_SOLUTIONS = true: the stored best schedules must stay consistent (asserted in checkBestSolutions);
# and checks that MIPs were solved, that math elimination never happens before effective elimination is possible
# for the last team, and that no mathematically eliminated team makes the playoffs (asserted in simulate).
# With Gurobi available, test/test_all.sh tests the same code with Gurobi.

using Tanking
using HiGHS

num_replications = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1

Tanking.MIP_OPTIMIZER[] = HiGHS.Optimizer
Tanking.eval(:(DEBUG = true))
Tanking.eval(:(CHECK_BEST_SOLUTIONS = true))
Tanking.set_mode(Tanking.STRICT)

failures = 0
for (step, seed) in [(1, 11), (16, 12), (31, 13)]
  t = @elapsed out = Tanking.simulate(Tanking.num_teams, Tanking.num_playoff_teams, 3, num_replications, Tanking.num_teams,
      0.71425, Tanking.breakpoint_list, Tanking.nba_odds_list, Tanking.nba_num_lottery, Tanking.true_strength,
      Tanking.STRICT, 2, [step], nothing, false; seed_per_step=seed, verbose=false)
  math_elim = out[7][step, end, 1]
  num_mips = out[9][step, 1]
  missing_case = out[17][step, 1]
  println("step $step ($(step-1) selfish teams): $(round(t, digits=1)) s; MIPs per replication = $num_mips; ",
      "teams mathematically eliminated at the end = $math_elim; num_missing_case = $missing_case")
  if !(num_mips > 0)
    println("FAIL: no MIPs were solved")
    global failures += 1
  end
  if !(0 < math_elim <= 30 - Tanking.num_playoff_teams)
    println("FAIL: number mathematically eliminated at the end of the season is $math_elim")
    global failures += 1
  end
end
println(failures == 0 ? "PASSED" : "FAILED", " (modes 2 and 3 agreed on every elimination check; best schedules consistent)")
exit(failures == 0 ? 0 : 1)
