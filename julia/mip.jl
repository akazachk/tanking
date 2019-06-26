########################################
# Mathematical elimination calculation #
########################################
# Aleksandr M. Kazachkov
# Shai Vardi
###
# Decide whether team k can still make the playoffs
#
# The MIP is used to calculate the best rank that team k can possibly still achieve
###

using JuMP, Cbc

function heuristicBestRank(k, t, schedule, in_stats, in_outcome,
    team_name_ind = 1, num_wins_ind = 2, games_left_ind = 3)
  ###
  # Heuristic for a ``good'' set of outcomes for team k for the remainder of the season
  # (remainder of the season = games t to num_games_total)
  ###
  num_games_total = size(schedule,1)
  num_teams = size(stats,1)

  ## Make copies of input data
  stats = copy(in_stats)
  outcome = copy(in_outcome)
  
  ## Calculate W
  thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome,
      team_name_ind, num_wins_ind, games_left_ind)
  W = stats[k, num_wins_ind]

  ## If a team's relative position to k is decided, it will win its remaining games
  heuristicHelper!(k, t, schedule, stats, outcome,
      team_name_ind, num_wins_ind, games_left_ind)

  ## Run rest of reason
  for game_ind = t:num_games_total
    # Skip the game if the outcome has already been decided
    if outcome[game_ind] > 0
      continue

    # Current teams playing
    i = schedule[game_ind,1] 
    j = schedule[game_ind,2]

    # Assign win to team with fewest wins
    # If same number of wins, assign it to the team with the _most_ games remaining
    if stats[i, num_wins_ind] < stats[j, num_wins_ind]
      winner = i
      loser = j
    else if stats[i, num_wins_ind] > stats[j, num_wins_ind]
      winner = j
      loser = i
    else if stats[i, games_left_ind] >= stats[j, games_left_ind]
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
      thisTeamWinsRemainingGames!(winner, game_ind, schedule, stats, outcome, 
          team_name_ind, num_wins_ind, games_left_ind)
    end
    if stats[loser, num_wins_ind]  + stats[loser, games_left_ind] <= W
      thisTeamWinsRemainingGames!(loser, game_ind, schedule, stats, outcome, 
          team_name_ind, num_wins_ind, games_left_ind)
    end
  end # iterate over remaining games

  ## Identify ranking of teams based on the outcomes
  sorted_teams = sortperm(stats[:,num_wins_ind], rev=true)
  rank_of_team = Array{Int}(undef, num_teams)
  rank_of_team[sorted_teams[1]] = 1
  for i = 2:num_teams
    last_team = sorted_teams[i-1]
    curr_team = sorted_teams[i]
    last_wins = w[last_team]
    curr_wins = w[curr_team]
    rank_of_team[last_team] = rank_of_team[curr_team]
    if last_wins > curr_wins
      rank_of_team[last_team] = rank_of_team[last_team] + 1
    end
  end

  return stats, outcome
end # heuristicBestRank

function thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome, 
    team_name_ind = 1, num_wins_ind = 2, games_left_ind = 3)
  ###
  # Set team k as winner of all its remaining games
  ###
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
      
    outcome[game_ind] = k
    stats[k, num_wins_ind] = stats[k, num_wins_ind] + 1
    for i in [k,j]
      stats[i,games_left_ind] = stats[i,games_left_ind] - 1
    end
  end # k wins remaining games
  return num_games_decided
end # thisTeamWinsRemainingGames

function heuristicHelper!(k, t, schedule, stats, outcome,
    team_name_ind = 1, num_wins_ind = 2, games_left_ind = 3)
  ###
  ## Teams with relative position to k decided will win their remaining games
  # * if i/j has more than W wins, that team wins
  # * else if i/j cannot have more than W wins, that team wins
  num_games_decided = 0
  game_ind = t
  W = stats[k, num_wins_ind]
  while game_ind <= num_games_total
    # Skip the game if the outcome has already been decided
    if outcome[game_ind] > 0
      continue

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
    else if (w_j > W)
      winner = j
    else if (W_i <= W)
      winner = i
    else if (W_j <= W)
      winner = j

    # If winner is determined, do updates for rest of season for that team, and reset game_ind
    if (winner > 0)
      num_games_decided += thisTeamWinsRemainingGames!(winner, game_ind, schedule, stats, outcome,
          team_name_ind, num_wins_ind, games_left_ind)
       game_ind = t-1
    end
    game_ind = game_ind + 1
  end
  return num_games_decided
end # heuristicHelper

function setupMIP(schedule, num_teams, num_team_games, num_games_total):
  ###
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
  model = Model(with_optimizer(Cbc.Optimizer))
  
  ## Set up variables and constraints
  @variable(model, w[1:num_teams]) # w_i = num wins of team i at end of season

  # x_{kt} = indicator that k wins game t
  game_ind_for_team = zeros(Int, num_teams)
  x_vars_for_team = Array(undef, num_teams, num_team_games)
  for t = 1:num_games_total
    i = schedule[t,1]
    j = schedule[t,2]
    for k in [i,j]
      game_ind_for_team[k] = game_ind_for_team[k] + 1
      x_vars_for_team[k,game_ind_for_team[k]] = 
          @variable(model, 
            base_name="x_{$k,$t}", 
            lower_bound = 0, upper_bound = 1,
            binary=true) # x_{kt} = indicator that k wins game t
    end
    xit = x_vars_for_team[k, game_ind_for_team[i]]
    xjt = x_vars_for_team[k, game_ind_for_team[j]]
    con = @constraint(model, xit + xjt = 1) # exactly one team wins game t
    set_name(con, "game$t")
  end # x variables

  @variable(model, y[1:num_teams] >= 0) # y_i = helper variable for linearizing max
  @variable(model, z[1:num_teams], lower_bound = 0, upper_bound = 1,
      binary=true) # z_i = indicator that team i has better rank than team k

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
  @objective(model, Min, 1 + sum(z[i] for i in 1:num_teams))

  return model
end # setupMIP

function resetMIP!(model, t, schedule, stats)
  ### 
  # Reset all outcomes after game t after fixing
  # This means reset all x_{it}, x_{jt}, and z variables
  ###
  num_games_total = size(schedule,1)
  num_teams = size(stats,1)

  fix(model[:W], 0)

  for game_ind = t:num_games_total
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    for k in [i,j]
      name = "x_{$k,$t}"
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
  else if j == k
    set_lower_bound(xjt, 1)
    set_upper_bound(xit, 0)
  end
  return
end # fixOutcome

function fixVariables!(model, k, t, W, schedule, stats)
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
end # fixVariables

function setIncumbent!(model, k, t, schedule, outcome)
  ### Set incumbent solution based on outcome of games for t,...,T ###
  num_games_total = size(schedule,1)
  num_teams = length(model[:z])

  W = 0
  w = zeros(Int, num_teams)
  x = -1 * ones(num_teams, num_games_total)
  y = zeros(Int, num_teams)
  z = zeros(Int, num_teams)

  ## Calculate number wins for each team in previous games
  for game_ind = 1:t-1
    i = schedule[game_ind,1]
    j = schedule[game_ind,2]
    xit = variable_by_name(model, "x_{$i,$t}")
    xjt = variable_by_name(model, "x_{$j,$t}")
    val_i = is_fixed(xit) ? fix_value(xit) : 0.5;
    val_j = is_fixed(xjt) ? fix_value(xjt) : 0.5;
    winner = 0
    if (val_i == 1)
      winner = i
    else if (val_j == 1)
      winner = j
    end
    w[winner] = w[winner] + 1
  end

  ## Set variables for remaining games and number of wins
  for game_ind = t:num_games_total
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

  return W, w, x, y, z
end # setIncumbent

function solveMIP!()
  return
end # solveMIP

function checkMIP(model)
  if termination_status(model) == MOI.OPTIMAL
    optimal_solution = value.(all_variables(model))
    optimal_objective = objective_value(model)
  elseif termination_status(model) == MOI.TIME_LIMIT && has_values(model)
    suboptimal_solution = value.(all_variables(model))
    suboptimal_objective = objective_value(model)
  else
    error("The model was not solved correctly.")
  end
  return
end # checkMIP
