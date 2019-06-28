########################################
# Mathematical elimination calculation #
########################################
# Aleksandr M. Kazachkov
# Shai Vardi
###
# Decide whether team k can still make the playoffs
#
# We calculate the best rank that team k can possibly still achieve
# A heuristic is provided, as well as a MIP formulation
# The heuristic can be used to provide a strong incumbent solution
###

using JuMP, MathOptFormat
#using Cbc
#using GLPK
using Gurobi

function heuristicBestRank(k, t, schedule, in_stats, in_outcome,
    num_wins_ind = 2, games_left_ind = 3)
  ###
  # Heuristic for a ``good'' set of outcomes for team k for the remainder of the season
  # (remainder of the season = games t to num_games_total)
  ###
  num_games_total = size(schedule,1)
  num_teams = size(in_stats,1)

  ## Make copies of input data
  stats = copy(in_stats)
  outcome = copy(in_outcome)
  
  ## Calculate W
  thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome, num_wins_ind, games_left_ind)
  W = stats[k, num_wins_ind]

  ## If a team's relative position to k is decided, it will win its remaining games
  heuristicHelper!(k, t, schedule, stats, outcome, num_wins_ind, games_left_ind)

  ## Run rest of reason
  for game_ind = t:num_games_total
    # Skip the game if the outcome has already been decided
    if outcome[game_ind] > 0
      continue
    end

    # Current teams playing
    i = schedule[game_ind,1] 
    j = schedule[game_ind,2]

    # Assign win to team with fewest wins
    # If same number of wins, assign it to the team with the _most_ games remaining
    if stats[i, num_wins_ind] < stats[j, num_wins_ind]
      winner = i
      loser = j
    elseif stats[i, num_wins_ind] > stats[j, num_wins_ind]
      winner = j
      loser = i
    elseif stats[i, games_left_ind] >= stats[j, games_left_ind]
      winner = i
      loser = j
    else
      winner = j
      loser = i
    end
    outcome[game_ind] = winner

    # Do updates
    stats[winner, num_wins_ind] = stats[winner, num_wins_ind] + 1
    stats[i,games_left_ind] = stats[i,games_left_ind] - 1
    stats[j,games_left_ind] = stats[j,games_left_ind] - 1

    # Check if more game outcomes can be fixed
    if stats[winner, num_wins_ind] > W
      thisTeamWinsRemainingGames!(winner, game_ind, schedule, stats, outcome, num_wins_ind, games_left_ind)
    end
    if stats[loser, num_wins_ind]  + stats[loser, games_left_ind] <= W
      thisTeamWinsRemainingGames!(loser, game_ind, schedule, stats, outcome, num_wins_ind, games_left_ind)
    end
  end # iterate over remaining games

  ## Identify ranking of teams based on the outcomes
  num_wins = stats[:,num_wins_ind]
  sorted_teams = sortperm(num_wins, rev=true)
  rank_of_team = Array{Int}(undef, num_teams)
  rank_of_team[sorted_teams[1]] = 1
  ct = 1
  for i = 2:num_teams
    last_team = sorted_teams[i-1]
    curr_team = sorted_teams[i]
    last_wins = num_wins[last_team]
    curr_wins = num_wins[curr_team]
    rank_of_team[curr_team] = rank_of_team[last_team]
    if last_wins > curr_wins
      rank_of_team[curr_team] += ct
      ct = 1
    else
      ct += 1
    end
  end

  return outcome, num_wins, rank_of_team
end # heuristicBestRank

function thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome, 
    num_wins_ind = 2, games_left_ind = 3)
  ###
  # Set team k as winner of all its remaining games
  ###
  num_games_total = length(outcome)
  num_games_decided = 0
  for game_ind = t:num_games_total
    if outcome[game_ind] > 0
      continue
    end
    if (schedule[game_ind, 1] != k) && (schedule[game_ind, 2] != k)
      continue
    end
    num_games_decided += 1

    if (schedule[game_ind, 1] != k)
      j = schedule[game_ind, 1]  
    else
      j = schedule[game_ind, 2]  
    end
      
    outcome[game_ind] = k
    stats[k, num_wins_ind] = stats[k, num_wins_ind] + 1
    for i in [k,j]
      stats[i,games_left_ind] = stats[i,games_left_ind] - 1
    end
  end # k wins remaining games
  return num_games_decided
end # thisTeamWinsRemainingGames

function heuristicHelper!(k, t, schedule, stats, outcome,
    num_wins_ind = 2, games_left_ind = 3)
  ###
  ## Teams with relative position to k decided will win their remaining games
  # * if i/j has more than W wins, that team wins
  # * elseif i/j cannot have more than W wins, that team wins
  num_games_total = length(outcome)
  num_games_decided = 0
  game_ind = t
  W = stats[k, num_wins_ind]
  while game_ind <= num_games_total
    # Skip the game if the outcome has already been decided
    if outcome[game_ind] > 0
      game_ind += 1
      continue
    end

    # Current teams playing
    i = schedule[game_ind,1] 
    j = schedule[game_ind,2]

    # Check whether winner can be set
    w_i = stats[i, num_wins_ind]
    g_i = stats[i, games_left_ind]
    W_i = w_i + g_i
    w_j = stats[i, num_wins_ind]
    g_j = stats[j, games_left_ind]
    W_j = w_j + g_j

    winner = 0
    if (w_i > W)
      winner = i
    elseif (w_j > W)
      winner = j
    elseif (W_i <= W)
      winner = i
    elseif (W_j <= W)
      winner = j
    end

    # If winner is determined, do updates for rest of season for that team, and reset game_ind
    if (winner > 0)
      num_games_decided += 
          thisTeamWinsRemainingGames!(winner, game_ind, schedule, stats, outcome, 
              num_wins_ind, games_left_ind)
       game_ind = t-1
    end
    game_ind += 1
  end
  return num_games_decided
end # heuristicHelper

function updateHeuristicBestRank!(winner, t, schedule, 
    best_outcomes, best_num_wins, best_rank)
  ###
  # Update set of best outcomes, num wins, ranks with outcome of game t
  ###
  num_teams = length(best_rank)
  loser = (winner == schedule[t,1]) ? schedule[t,2] : schedule[t,1]

  for i = 1:num_teams
    # Skip the uninitialized teams, those for which the outcome matches, and elimiinated teams
    if best_rank[i] <= 0 || best_outcomes[i,t] == winner
      continue
    end
    
    # For the remaining teams, the outcome does not match
    # This means the winner has one extra win at the end of the season,
    # and this may change the best possible rank the other teams can achieve
    best_outcomes[i,t] = winner
    best_num_wins[i, winner] += 1
    best_num_wins[i, loser] -= 1
    if (i != winner && i != loser)
      if (best_num_wins[i, winner] - 1 == best_num_wins[i, i])
        best_rank[i] = best_rank[i] + 1 # winner moves out of team i's equivalence class
      end
      if (best_num_wins[i, loser] == best_num_wins[i, i])
        best_rank[i] = best_rank[i] - 1 # loser moved into team i's equivalence class
      end
    else
      # Count number of teams with strictly more wins; rank is 1 greater than this value
      val = count(j->(j > best_num_wins[i,i]), best_num_wins[i,:])
      best_rank[i] = val + 1
    end
  end

  return
end # updateHeuristicBestRank

function setupMIP(schedule, num_teams, num_playoff_teams, num_team_games, num_games_total)
  ###
  # Arguments:
  #   schedule: Matrix{Int}(num_games_total, 2)
  #   num_teams, num_playoff_teams, num_team_games, num_games_total: scalars
  #
  # Return the following MIP model
  #
  # Constants:
  #   k: index of the team we are solving the min rank problem for
  #   M: total number of games each team plays
  #   G_t: set of two teams that are playing in game t (for each t \in [T])
  #
  # Variables:
  #   W: initial num wins of team k + num games remaining for k (set to 0 now)
  #   w_i: continuous; num wins of team i at end of season
  #   x_{it}: binary; whether team i wins game t
  #   y_i: continuous; represents max{0, w_i - W}
  #   z_i: binary; whether team i is ranked above team k at end
  #
  # Objective:
  #   min 1 + \sum_i z_i
  #
  # Constraints:
  #   \sum_{i \in G_t} x_{it} = 1                         (for all t)
  #   w_i = \sum_t x_{it}                                 (for all i)
  #   z_i \ge (w_i - W) / M                               (for all i)
  #   y_i \ge z_i                                         (for all i)
  #   y_i \ge w_i - W                                     (for all i)
  #
  # Additional constraints with respect to k (not added here):
  #   w_k = W
  #   x_{kt} = 1 for all t in which team k plays
  #   y_k = 0
  #   z_k = 0
  #
  # Binaries:
  #   x, z \in \{0,1\}
  #
  # Bounds:
  #   y_i \ge 0                                           (for all i)
  ###
  #model = Model(with_optimizer(Cbc.Optimizer, logLevel=0)) # about five times slower than Gurobi (or worse)
  #model = Model(with_optimizer(GLPK.Optimizer))
  model = Model(with_optimizer(Gurobi.Optimizer, BestObjStop=num_playoff_teams, BestBdStop=num_playoff_teams, TimeLimit=10, OutputFlag=0))
  
  ## Set up variables and constraints
  @variable(model, w[1:num_teams]) # w_i = num wins of team i at end of season

  # x_{kt} = indicator that k wins game t
  @variable(model, 0 <= x_vars_for_team[1:num_teams,1:num_team_games] <= 1, Int)
  game_ind_for_team = zeros(Int, num_teams)
  for t = 1:num_games_total
    i = schedule[t,1]
    j = schedule[t,2]
    for k in [i,j]
      game_ind_for_team[k] = game_ind_for_team[k] + 1
    end
    xit = x_vars_for_team[i, game_ind_for_team[i]]
    xjt = x_vars_for_team[j, game_ind_for_team[j]]
    set_name(xit, "x_{$i,$t}")
    set_name(xjt, "x_{$j,$t}")
    con = @constraint(model, xit + xjt == 1) # exactly one team wins game t
    set_name(con, "game$t")
  end # x variables

  @variable(model, y[1:num_teams] >= 0) # y_i = helper variable for linearizing max
  @variable(model, z[1:num_teams], lower_bound = 0, upper_bound = 1,
      integer=true) # z_i = indicator that team i has better rank than team k

  # Set team-wise constraints
  @variable(model, W == 0)
  for i = 1:num_teams
    con = @constraint(model, w[i] - sum(x_vars_for_team[i,:]) == 0) # num wins of team i at season end
    set_name(con, "w$i")
    con = @constraint(model, num_team_games * z[i] >= w[i] - W)
    set_name(con, "LBz$i")
    con = @constraint(model, y[i] - z[i] >= 0)
    set_name(con, "UBz$i")
    con = @constraint(model, y[i] >= w[i] - W)
    set_name(con, "LBy$i")
  end

  ## Add objective
  @variable(model, rank)
  @constraint(model, rk, rank == sum(z) + 1)

  ## Make it an optimization problem or feasibility problem
  @objective(model, Min, rank)
  #@constraint(model, cutoff, rank <= num_playoff_teams)

  return model
end # setupMIP

function resetMIP!(model, t, schedule, stats)
  ### 
  # Reset all outcomes after game t after fixing
  # This means reset W, x_{it}, x_{jt}, and z variables
  ###
  num_games_total = size(schedule,1)
  num_teams = size(stats,1)

  fix(model[:W], 0)

  for game_ind = t:num_games_total
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    for k in [i,j]
      name = "x_{$k,$game_ind}"
      xkt = variable_by_name(model, name)
      set_lower_bound(xkt, 0)
      set_upper_bound(xkt, 1)
    end
  end

  for i = 1:num_teams
    z = variable_by_name(model, "z[$i]")
    set_lower_bound(z, 0)
    set_upper_bound(z, 1)
  end

  return
end # resetMIP

function fixOutcome!(model, k, t, schedule)
  ###
  # Winner of game t is now determined, so set set x_it and x_jt appropriately
  ###
  i = schedule[t,1]
  j = schedule[t,2]
  xit = variable_by_name(model, "x_{$i,$t}")
  xjt = variable_by_name(model, "x_{$j,$t}")
  if i == k
    set_lower_bound(xit, 1)
    set_upper_bound(xjt, 0)
  elseif j == k
    set_lower_bound(xjt, 1)
    set_upper_bound(xit, 0)
  end

  return
end # fixOutcome

function fixVariables!(model, k, t, W, schedule)
  ###
  # Fix variables with respect to team k
  #
  # This means adjusting 
  #   * W (i.e., changing the right-hand side of the appropriate constraints)
  #   * unfixing old variables that were fixed (change their bounds back to [0,1])
  #   * fixing new variables as needed
  ###
  num_games_total = size(schedule,1)
  
  fix(model[:W], W)

  for game_ind = t:num_games_total
    fixOutcome!(model, k, game_ind, schedule)
  end
  z = variable_by_name(model, "z[$k]")
  set_upper_bound(z, 0)

  return
end # fixVariables

function setIncumbent!(model, k, t, schedule, outcome)
  ### 
  # Set incumbent solution based on given outcome of games for 1,...,T
  # It is assumed that the outcomes of games 1,...,t-1 match the fixed
  # value of the game outcomes in model for those games
  ###
  num_games_total = size(schedule,1)
  num_teams = length(model[:z])

  W = 0
  w = zeros(Int, num_teams)
  x = -1 * ones(num_teams, num_games_total)
  y = zeros(Int, num_teams)
  z = zeros(Int, num_teams)

  ## Calculate number wins for each team in previous games,
  ## and set variables for remaining games and number of wins
  for game_ind = 1:num_games_total
    winner = outcome[game_ind]
    if winner <= 0
      continue
    end
    w[winner] = w[winner] + 1
    
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    if i == winner
      x[i,game_ind] = 1
      x[j,game_ind] = 0
    else
      x[j,game_ind] = 1
      x[i,game_ind] = 0
    end
  end

  W = w[k]

  ## Set y and z based on the final outcomes
  for i = 1:num_teams
    if w[i] > W
      y[i] = w[i] - W
      z[i] = 1
    else
      y[i] = 0
      z[i] = 0
    end
  end

  ## Set start values
  # W
  set_start_value(model[:W], W)
  # x
  for game_ind = t:num_games_total
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    xit = variable_by_name(model, "x_{$i,$game_ind}")
    xjt = variable_by_name(model, "x_{$j,$game_ind}")
    set_start_value(xit, x[i,game_ind])
    set_start_value(xjt, x[j,game_ind])
  end
  for i = 1:num_teams
    # w
    set_start_value(model[:w][i], w[i])
    # y
    set_start_value(model[:y][i], y[i])
    # z
    set_start_value(model[:z][i], z[i])
  end

  return W, w, x, y, z
end # setIncumbent

function solveMIP!(model, num_playoff_teams)
  optimize!(model)
  status = termination_status(model)
  if status == MOI.OPTIMAL
    return true
  elseif status == MOI.OBJECTIVE_LIMIT
    #if objective_value(model) <= num_playoff_teams
    if value.(model[:rank]) <= num_playoff_teams
      return true
    else
      return false
    end
  elseif status == MOI.TIME_LIMIT
    ## Save the hard LP
    lp_file = MathOptFormat.LP.Model()
    MOI.copy_to(lp_file, backend(model))
    MOI.write_to_file(lp_file, "hard.lp")

    if has_values(model)
      #if objective_value(model) <= num_playoff_teams
      if value.(model[:rank]) <= num_playoff_teams
        return true
      else
        return false
      end
    else
      return false
    end
  elseif status == MOI.INFEASIBLE || status == MOI.INFEASIBLE_OR_UNBOUNDED
    return false
  else
    error("The model was not solved correctly. Exiting with status $status.")
  end
  return false
end # solveMIP

function updateUsingMIPSolution!(model, k, t, schedule,
    best_outcomes, best_num_wins, best_rank)
  status = termination_status(model)
  good_status = (status == MOI.OPTIMAL) || (status == MOI.OBJECTIVE_LIMIT)
  if !good_status || !has_values(model)
    error("The model was not solved correctly. Exiting with status $status.")
  end

  num_games_total = size(schedule, 1)
  num_teams = size(best_rank, 1)

  for game_ind = t:num_games_total
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    xit = variable_by_name(model, "x_{$i,$game_ind}")
    xjt = variable_by_name(model, "x_{$j,$game_ind}")
    winner = isVal(value.(xit), 1) ? i : j

    best_outcomes[k, game_ind] = winner
  end

  w = model[:w]
  for i = 1:num_teams
    best_num_wins[k, i] = Int(round(value.(w[i])))
  end

  #best_rank[k] = Int(round(objective_value(model)))
  best_rank[k] = Int(round(value.(model[:rank])))
end # updateUsingMIPSolution

function updateOthersUsingBestSolution!(k, t, schedule,
    best_outcomes, best_num_wins, best_rank)
  ###
  # Check whether other teams best schedule can be updated
  ###
  num_games_total = length(schedule)
  num_teams = length(best_rank)

  num_wins = best_num_wins[k,:]
  sorted_teams = sortperm(num_wins, rev=true)
  rank_of_team = Array{Int}(undef, num_teams)
  rank_of_team[sorted_teams[1]] = 1
  ct = 1
  for i = 2:num_teams
    last_team = sorted_teams[i-1]
    curr_team = sorted_teams[i]
    last_wins = num_wins[last_team]
    curr_wins = num_wins[curr_team]
    rank_of_team[curr_team] = rank_of_team[last_team]
    if last_wins > curr_wins
      rank_of_team[curr_team] += ct
      ct = 1
    else
      ct += 1
    end
  end

  for i = 1:num_teams
    if best_rank[i] == 0 || best_rank[i] > rank_of_team[i]
      best_outcomes[i, :] = best_outcomes[k, :]
      best_num_wins[i, :] = best_num_wins[k, :]
      best_rank[i] = rank_of_team[i]
    end
  end
end # updateOthersUsingBestSolution
