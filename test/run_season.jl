using Tanking # "using" adds symbols into namespace
import Gurobi  # for Gurobi.Env()

global GRB_ENV = Gurobi.Env()

#import Random # "import" keeps namespace clean
#Random.seed!(628) # for reproducibility; NOTE: should be set once somewhere before this function is run

function run_season(;true_strength=nothing)

  num_teams = 30
  num_playoff_teams = Int(2^ceil(log(2, num_teams / 2)))
  num_rounds = 3
  num_replications = 1
  num_steps = 1
  gamma = 0.71425
  breakpoint_list = [1.]
  nba_odds_old = [.250, .199, .156, .119, .088, .063, .043, .028, .017, .011, .008, .007, .006, .005]
  nba_odds_old = nba_odds_old[14:-1:1]
  nba_odds_new = [.140, .140, .140, .125, .105, .090, .075, .060, .045, .030, .020, .015, .010, .005]
  nba_odds_new = nba_odds_new[14:-1:1]
  nba_odds_flat = [1. / 14. for i in 1:14]
  nba_odds_list = [nba_odds_new, nba_odds_old, nba_odds_flat]
  nba_num_lottery = [4, 3, 4]
  if isnothing(true_strength)
    true_strength = num_teams:-1:1
  end
  math_elim_mode = 0
  mode = Tanking.STRICT
  mode_list_name = ["BT"]
  selected_steps = 1

  print("\n## Running one season. ##\nOutput is win_pct[step_ind, team_ind, stat], where step_ind is an integer corresponding to a repeat of the simulation with a different ratio of teams tanking, team_ind takes a value from 1 to num_teams, and stat is either avg_stat = 1, stddev_stat = 2, min_stat = 3, max_stat = 4.\n\n")
  win_pct = Tanking.simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, nba_num_lottery, true_strength, mode, math_elim_mode, selected_steps, GRB_ENV, true)

  ## Stats we keep
  avg_stat    = 1
  stddev_stat = 2
  min_stat    = 3
  max_stat    = 4
  num_stats   = 4
  prefix      = ["avg_", "stddev_", "min_", "max_"]

  display(win_pct[1, :, avg_stat])
end
