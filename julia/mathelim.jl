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

"""
heuristicBestRank: Heuristic for a ``good'' set of outcomes for team k for the remainder of the season

In this heuristic, we assume team k wins its remaining games
(remainder of the season = games t to num_games_total)

Say this means team k ends up with W wins
We first check whether there is already a schedule we can copy
  -> if there is a schedule in which the last playoff team has W or fewer wins
Any team that can never or will always have W wins (or more) will win their remaining games

Returns
  * set of outcomes
  * head-to-head records
  * number of wins each team has
  * rank of each team
"""
function heuristicBestRank(k, t, schedule, in_stats, in_outcome, in_h2h,
    best_outcomes, best_h2h, best_num_wins, best_rank, num_playoff_teams,
    num_wins_ind = 2, games_left_ind = 3)
  num_games_total = size(schedule,1)
  num_teams = size(in_stats,1)

  ## Check we are not immediately done
  if best_rank[k] > 0 && best_rank[k] <= num_playoff_teams
    outcome = best_outcomes[k,:]
    h2h = best_h2h[k,:,:]
    num_wins = best_num_wins[k,:]
    rank_of_team = rankTeamsFromWinTotals(num_wins)
    return outcome, h2h, num_wins, rank_of_team
  end
  
  ## Make copies of input data
  # TODO is calling copy needed?
  stats = copy(in_stats)
  outcome = copy(in_outcome)
  h2h = copy(in_h2h)
  
  ## Calculate W
  thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome, h2h, num_wins_ind, games_left_ind)
  W = stats[k, num_wins_ind]

  ## Find whether we can copy an existing schedule
  for i = 1:num_teams
    if i == k || best_rank[i] <= 0
      continue
    end
    sorted_teams = sortperm(best_num_wins[i,:], rev=true)
    W_i = best_num_wins[i,sorted_teams[num_playoff_teams]]
    if W_i <= W
      # Found a schedule to use
      # Set outcome of all remaining games in which k does not play based on best_outcomes[i,:]
      for game_ind = t:num_games_total
        # Skip the game if the outcome has already been decided
        if outcome[game_ind] > 0
          continue
        end

        # Current teams playing
        i = schedule[game_ind,1] 
        j = schedule[game_ind,2]

        # Assign win
        outcome[game_ind] = best_outcomes[i,game_ind]
        winner = outcome[game_ind]
        loser = (winner == i) ? i : j

        # Do updates
        stats[winner, num_wins_ind] += 1
        stats[i,games_left_ind] -= 1
        stats[j,games_left_ind] -= 1
        h2h[winner,loser] += 1
      end # iterate over remaining games
      num_wins = stats[:,num_wins_ind]
      rank_of_team = rankTeamsFromWinTotals(num_wins)
      return outcome, h2h, num_wins, rank_of_team
    end
  end

  ## If a team's relative position to k is decided, it will win its remaining games
  heuristicHelper!(k, t, schedule, stats, outcome, h2h, num_wins_ind, games_left_ind)

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
    stats[winner, num_wins_ind] += 1
    stats[i,games_left_ind] -= 1
    stats[j,games_left_ind] -= 1
    h2h[winner,loser] += 1

    # Check if more game outcomes can be fixed
    if stats[winner, num_wins_ind] > W
      thisTeamWinsRemainingGames!(winner, game_ind, schedule, 
          stats, outcome, h2h, num_wins_ind, games_left_ind)
    end
    if stats[loser, num_wins_ind] + stats[loser, games_left_ind] <= W
      thisTeamWinsRemainingGames!(loser, game_ind, schedule, 
          stats, outcome, h2h, num_wins_ind, games_left_ind)
    end
  end # iterate over remaining games

  ## Identify ranking of teams based on the outcomes
  num_wins = stats[:,num_wins_ind]
  rank_of_team = rankTeamsFromWinTotals(num_wins)

  return outcome, h2h, num_wins, rank_of_team
end # heuristicBestRank

"""
thisTeamWinsRemainingGames!: Set team k as winner of all its remaining games

Updates stats, outcome, and h2h
"""
function thisTeamWinsRemainingGames!(k, t, schedule, stats, outcome, h2h,
    num_wins_ind = 2, games_left_ind = 3)
  num_games_total = length(outcome)
  num_games_decided = 0
  for game_ind = t:num_games_total
    if outcome[game_ind] > 0
      continue # outcome is already decided and accounted for
    end
    if (schedule[game_ind, 1] != k) && (schedule[game_ind, 2] != k)
      continue # only look at games in which team k plays
    end
    num_games_decided += 1

    if (schedule[game_ind, 1] != k)
      j = schedule[game_ind, 1]  
    else
      j = schedule[game_ind, 2]  
    end
      
    outcome[game_ind] = k
    stats[k, num_wins_ind] += 1
    h2h[k,j] += 1
    for i in [k,j]
      stats[i,games_left_ind] -= 1
    end
  end # k wins remaining games
  return num_games_decided
end # thisTeamWinsRemainingGames

"""
heuristicHelper!: Teams with relative position to k decided will win their remaining games

    if i/j has more than W wins, that team wins
    elseif i/j cannot have more than W wins, that team wins

Updates stats, outcome, h2h
"""
function heuristicHelper!(k, t, schedule, stats, outcome, h2h,
    num_wins_ind = 2, games_left_ind = 3)
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
          thisTeamWinsRemainingGames!(winner, game_ind, schedule, 
              stats, outcome, h2h, num_wins_ind, games_left_ind)
       game_ind = t-1
    end
    game_ind += 1
  end
  return num_games_decided
end # heuristicHelper

"""
updateHeuristicBestRank!: Update set of best outcomes, num wins, ranks with outcome of game t

We assume that h2h has already been updated with the win
"""
function updateHeuristicBestRank!(winner, t, schedule, h2h,
    best_outcomes, best_h2h, best_num_wins, best_rank)
  num_teams = length(best_rank)
  num_games_total = size(schedule,1)
  loser = (winner == schedule[t,1]) ? schedule[t,2] : schedule[t,1]

  for i = 1:num_teams
    ## Skip the uninitialized teams, those for which the outcome matches, and eliminated teams
    if best_rank[i] <= 0 || best_outcomes[i,t] == winner
      continue
    end

    ## Check whether we can reallocate this win
    ## Our previous best schedule had the loser winning, instead of the winner
    ## If the winner has another game it was supposed to win in the best schedule,
    ## then we can "swap" that future win with a win now
    ## Note the -1 because we assume that h2h has already been updated
    num_future_wins_by_winner = best_h2h[i,winner,loser] - (h2h[winner,loser] - 1)
    if num_future_wins_by_winner > 0
      # Find the next game that winner was supposed to win
      game_ind = t+1
      while game_ind <= num_games_total
        if best_outcomes[i,game_ind] == winner
          break
        end
        game_ind += 1
      end

      # Set winner as winner of current game, and loser as winner of game_ind
      if game_ind <= num_games_total
        best_outcomes[i,t] = winner
        best_outcomes[i,game_ind] = loser
      end
      continue # continue iterating through the teams
    end
    
    ## For the remaining teams, the outcome does not match
    ## This means the winner has one best win at the end of the season,
    ## and this may change the best possible rank the other teams can achieve
    best_outcomes[i,t] = winner
    best_h2h[i,winner,loser] += 1
    best_h2h[i,loser,winner] -= 1
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

"""
setupMIP

`math_elim_mode`: <= 0: do not use math elim; 1: do not use MIP for math elim; 2-3: team-wise formulation; 4-5: cutoff formulation
"""
function setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
  if math_elim_mode < 1
    return 0
  elseif math_elim_mode in [2,3]
    return setupMIPByTeam(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
  elseif math_elim_mode in [4,5]
    return setupMIPByCutoff(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
  end
end # setupMIP

"""
SetupMIPByTeam

Parameters
---
   `schedule`: Matrix{Int}(num_games_total, 2)
   `h2h_left`: Matrix{Int}{num_teams,num_teams}
   `num_teams, `num_playoff_teams`, `num_team_games`, `num_games_total`: scalars
   `math_elim_mode`: see main.jl or simulate.jl

Returns
---
Return the following MIP model

```
 Constants:
   `k`: index of the team we are solving the min rank problem for
   `M`: total number of games each team plays
   `G_t`: set of two teams that are playing in game `t` (for each `t \\in [T]`)

 Variables:
   `W`: initial num wins of team `k` + num games remaining for `k` (set to 0 now)
   `w_i`: continuous; num wins of team `i` at end of season
   `math_elim_mode` == 2: x_{it}: binary; whether team `i` wins game `t`
   `math_elim_mode` == 3: x_{ij}: general integer; number of wins team `i` has over team `j`
   `y_i`: continuous; represents `max{0, w_i - W}`
   `z_i`: binary; whether team `i` is ranked above team `k` at end

 Objective:
   min rank = 1 + \\sum_i z_i

 Constraints:
   math_elim_mode == 2: \\sum_{i \\in G_t} x_{it} = 1    (for all t)
   math_elim_mode == 3: x_{ij} + x_{ji} = g_{ij}       (for all i,j)
   w_i = \\sum x_{i,:}                                  (for all i)
   z_i \\ge (w_i - W) / M                               (for all i)
   y_i \\ge z_i                                         (for all i)
   y_i \\ge w_i - W                                     (for all i)

 Binaries:
   math_elim_mode == 2: x \\in \\{0,1\\}
   z \\in \\{0,1\\}

 Integers:
   math_elim_mode == 3: x \\in [0,g_{ij}]

 Bounds:
   x \\ge 0
   math_elim_mode == 3: x_{i,i} = 0
   y_i \\ge 0                                           (for all i)
```
"""
function setupMIPByTeam(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
  if math_elim_mode != 2 && math_elim_mode != 3
    return 0
  end
  
  #model = Model(with_optimizer(Cbc.Optimizer, logLevel=0)) # about five times slower than Gurobi (or worse)
  #model = Model(with_optimizer(GLPK.Optimizer))
  model = Model(with_optimizer(Gurobi.Optimizer, BestObjStop=num_playoff_teams+1e-3, BestBdStop=num_playoff_teams+1e-3, TimeLimit=10, OutputFlag=0))
  
  ## Set up variables and constraints
  @variable(model, w[1:num_teams]) # w_i = num wins of team i at end of season

  if math_elim_mode == 2
    # x_{it} = indicator that team i wins game t
    @variable(model, 0 <= x[1:num_teams,1:num_team_games] <= 1, Int)
    game_ind_for_team = zeros(Int, num_teams)
    for t = 1:num_games_total
      i = schedule[t,1]
      j = schedule[t,2]
      for k in [i,j]
        game_ind_for_team[k] = game_ind_for_team[k] + 1
      end
      xit = x[i, game_ind_for_team[i]]
      xjt = x[j, game_ind_for_team[j]]
      #str = "x[$i,$t]"
      #xit = @variable(model, base_name=str, lower_bound=0, upper_bound=1, integer=true)
      #str = "x[$j,$t]"
      #xjt = @variable(model, base_name=str, lower_bound=0, upper_bound=1, integer=true)
      #xit = @variable(model, lower_bound=0, upper_bound=1, integer=true)
      #xjt = @variable(model, lower_bound=0, upper_bound=1, integer=true)
      set_name(xit, "x[$i,$t]")
      set_name(xjt, "x[$j,$t]")
      con = @constraint(model, xit + xjt == 1) # exactly one team wins game t
      set_name(con, "game$t")
    end
  elseif math_elim_mode == 3
    # x_{ij} = number games that team i wins over team j
    @variable(model, x[1:num_teams,1:num_teams] >= 0, Int)
    for i = 1:num_teams
      set_upper_bound(x[i,i], 0)
      for j = 1:num_teams
        if i == j
          continue
        end
        con = @constraint(model, x[i,j] + x[j,i] == h2h_left[i,j]) # total wins in series matches total games played
        set_name(con, "series_$(i)_$(j)")
        set_upper_bound(x[i,j], h2h_left[i,j]) # redundant, but could be useful for some solvers
        set_upper_bound(x[j,i], h2h_left[j,i]) # redundant, but could be useful for some solvers
      end
    end
  end # x variables

  @variable(model, y[1:num_teams] >= 0) # y_i = helper variable for linearizing max
  @variable(model, z[1:num_teams], lower_bound = 0, upper_bound = 1,
      integer=true) # z_i = indicator that team i has better rank than team k

  ## Set team-wise constraints
  @variable(model, W == 0)
  for i = 1:num_teams
    con = @constraint(model, w[i] - sum(x[i,:]) == 0) # num wins of team i at season end
    set_name(con, "w$i")
    con = @constraint(model, num_team_games * z[i] >= w[i] - W)
    set_name(con, "LBz$i")
    con = @constraint(model, y[i] - z[i] >= 0)
    set_name(con, "UBz$i")
    con = @constraint(model, y[i] >= w[i] - W)
    set_name(con, "LBy$i")
  end

  ## Add the objective
  @variable(model, rank)
  @constraint(model, rk, rank == sum(z) + 1)

  ## We can make it an optimization problem (adding the min rank objective) or a feasibility problem (does there exist a solution in which rank <= num_playoff_teams)
  @objective(model, Min, rank)
  #@constraint(model, cutoff, rank <= num_playoff_teams)

  return model
end # setupMIPByTeam

"""
setupMIPByCutoff

Parameters
---
* `schedule`: Matrix{Int}(num_games_total, 2)
* `h2h_left`: Matrix{Int}{num_teams,num_teams}
* `num_teams`, `num_playoff_teams`, `num_team_games`, `num_games_total`: scalars
* `math_elim_mode`: see main.jl or simulate.jl

Returns
---
Return the following MIP model

```
  Constants:
    k: index of the team we are solving the min rank problem for
    M: total number of games each team plays
    G_t: set of two teams that are playing in game t (for each t \\in [T])
    n^*: number of teams in playoffs
 
  Variables:
    W: number of wins by last team that makes the playoffs
    w_i: number of wins by team i
    math_elim_mode == 4: x_{it}: binary; whether team i wins game t
    math_elim_mode == 5: x_{ij}: general integer; number of wins team i has over team j
    alpha_i: binary; 0 if W >= num wins of team i (i.e., will = 1 for first n^* teams)
 
  Objective:
    min W
 
  Constraints:
    math_elim_mode == 4: \\sum_{i \\in G_t} x_{it} = 1    (for all t)
    math_elim_mode == 5: x_{ij} + x_{ji} = g_{ij}       (for all i,j)
    w_i = \\sum x_{i,:}
    W \\ge \\sum x_{i,:} - M \\alpha_i                     (for all i)
    \\sum_i \\alpha_i = n^*                               (for all i)
 
  Binaries:
    math_elim_mode == 4: x \\in \\{0,1\\}
    alpha \\in {0,1}^n
 
  Integers:
    math_elim_mode == 5: x_{ij} \\in [0,g_{ij}]
 
  Bounds:
    math_elim_mode == 5: x_{ii} = 0
```
"""
function setupMIPByCutoff(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
  if math_elim_mode != 4 && math_elim_mode != 5
    return 0
  end
  
  #model = Model(with_optimizer(Cbc.Optimizer, logLevel=0)) # about five times slower than Gurobi (or worse)
  #model = Model(with_optimizer(GLPK.Optimizer))
  #model = Model(with_optimizer(Gurobi.Optimizer, BestObjStop=num_playoff_teams, BestBdStop=num_playoff_teams, TimeLimit=10, OutputFlag=0))
  model = Model(with_optimizer(Gurobi.Optimizer, TimeLimit=10, OutputFlag=0))
  
  ## Set up variables and constraints
  @variable(model, W >= 0)
  @variable(model, w[1:num_teams]) # w_i = num wins of team i at end of season
  @variable(model, alpha[1:num_teams], lower_bound = 0, upper_bound = 1,
      integer=true) # alpha_i = indicator that team i has better rank than n^*
  con = @constraint(model, alpha_bd, sum(alpha) == num_playoff_teams)

  if math_elim_mode == 4
    # x_{it} = indicator that team i wins game t
    @variable(model, 0 <= x[1:num_teams,1:num_team_games] <= 1, Int)
    game_ind_for_team = zeros(Int, num_teams)
    for t = 1:num_games_total
      i = schedule[t,1]
      j = schedule[t,2]
      for k in [i,j]
        game_ind_for_team[k] = game_ind_for_team[k] + 1
      end
      xit = x[i, game_ind_for_team[i]]
      xjt = x[j, game_ind_for_team[j]]
      set_name(xit, "x[$i,$t]")
      set_name(xjt, "x[$j,$t]")
      con = @constraint(model, xit + xjt == 1) # exactly one team wins game t
      set_name(con, "game$t")
    end
  elseif math_elim_mode == 5
    # x_{ij} = number games that team i wins over team j
    @variable(model, x[1:num_teams,1:num_teams] >= 0, Int)
    for i = 1:num_teams
      set_upper_bound(x[i,i], 0)
      for j = 1:num_teams
        if i == j
          continue
        end
        con = @constraint(model, x[i,j] + x[j,i] == h2h_left[i,j]) # total wins in series matches total games played
        set_name(con, "series_$(i)_$(j)")
        set_upper_bound(x[i,j], h2h_left[i,j]) # redundant, but could be useful for some solvers
        set_upper_bound(x[j,i], h2h_left[j,i]) # redundant, but could be useful for some solvers
      end
    end
  end # x variables

  ## Set team-wise constraints
  for i = 1:num_teams
    con = @constraint(model, w[i] - sum(x[i,:]) == 0) # num wins of team i at season end
    set_name(con, "w$i")
    con = @constraint(model, w[i] - num_team_games * alpha[i] <= W) # bound on num wins of last playoff team
    set_name(con, "W$i")
  end

  ## Add the objective
  @objective(model, Min, W)

  return model
end # setupMIPByCutoff

function resetMIP!(model, t, schedule, h2h, stats, math_elim_mode)
  ### 
  # Reset all outcomes after game t after fixing
  # This means reset W, x_{it}, x_{jt}, and z variables
  # (or x_{ij} and x_{ji} variables in general integer case)
  ###
  num_games_total = size(schedule,1)
  num_teams = size(stats,1)

  if math_elim_mode in [2,4]
    for game_ind = t:num_games_total
      i = schedule[game_ind,1]
      j = schedule[game_ind,2]
      for k in [i,j]
        name = "x[$k,$game_ind]"
        xkt = variable_by_name(model, name)
        set_lower_bound(xkt, 0)
        set_upper_bound(xkt, 1)
      end
    end
  elseif math_elim_mode in [3,5]
    for i = 1:num_teams
      for j = 1:num_teams
        if i == j
          continue
        end
        set_lower_bound(variable_by_name(model, "x[$i,$j]"), h2h[i,j])
        set_lower_bound(variable_by_name(model, "x[$j,$i]"), h2h[j,i])
      end
    end
  end

  if math_elim_mode in [2,3]
    fix(model[:W], 0)
    for i = 1:num_teams
      z = variable_by_name(model, "z[$i]")
      set_lower_bound(z, 0)
      set_upper_bound(z, 1)
    end
  end

  return
end # resetMIP

function fixOutcome!(model, k, t, schedule, math_elim_mode)
  ###
  # Winner of game t is now determined, so set set x_it and x_jt appropriately
  # (or x_ij, x_jt, for the non-binary formulation)
  ###
  if math_elim_mode <= 1 # only update when non-heuristic methods are used
    return
  end
  i = schedule[t,1]
  j = schedule[t,2]
  if (k != i) && (k != j)
    return
  end

  if math_elim_mode in [2,4]
    xit = variable_by_name(model, "x[$i,$t]")
    xjt = variable_by_name(model, "x[$j,$t]")
    if i == k
      set_lower_bound(xit, 1)
      set_upper_bound(xjt, 0)
    elseif j == k
      set_lower_bound(xjt, 1)
      set_upper_bound(xit, 0)
    end
  end # binary case

  if math_elim_mode in [3,5]
    xij = variable_by_name(model, "x[$i,$j]")
    xji = variable_by_name(model, "x[$j,$i]")
    lb_ij = lower_bound(xij)
    lb_ji = lower_bound(xji)
    if i == k
      set_lower_bound(xij, lb_ij + 1)
    elseif j == k
      set_lower_bound(xji, lb_ji + 1)
    end
  end # general integer case

  return
end # fixOutcome

"""
fixVariables!
Fix variables with respect to team k

We do not do this for the cutoff-based formulation

This means adjusting 
  * W (i.e., changing the right-hand side of the appropriate constraints)
  * fixing new variables as needed
"""
function fixVariables!(model, k, t, W, schedule, math_elim_mode)
  if math_elim_mode < 1 || math_elim_mode in [4,5]
    return
  end

  num_games_total = size(schedule,1)
  
  for game_ind = t:num_games_total
    fixOutcome!(model, k, game_ind, schedule, math_elim_mode)
  end

  if math_elim_mode in [2,3]
    fix(model[:W], W)
    z = variable_by_name(model, "z[$k]")
    set_upper_bound(z, 0)
  end

  return
end # fixVariables

function setIncumbent!(model, k, t, schedule, outcome, math_elim_mode)
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
    xit = variable_by_name(model, "x[$i,$game_ind]")
    xjt = variable_by_name(model, "x[$j,$game_ind]")
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

function solveMIP!(model)
  ###
  # Return true if has values that can be used
  ###
  optimize!(model)
  status = termination_status(model)
  good_status = status in [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.TIME_LIMIT]
  other_status = status in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
  if !(good_status || other_status)
    error("The model was not solved correctly. Exiting with status $status.")
  end
  return good_status && has_values(model)
end # solveMIP

function checkMIP(model, cutoff, math_elim_mode)
  ###
  # Check whether MIP says the team is eliminated
  # If there is a time out, it may give the incorrect answer
  ###
  status = termination_status(model)
  if status == MOI.OPTIMAL
    return true
  elseif status in [MOI.OBJECTIVE_LIMIT, MOI.TIME_LIMIT]
    if status == MOI.TIME_LIMIT
      ## Save the hard LP
      lp_file = MathOptFormat.LP.Model()
      MOI.copy_to(lp_file, backend(model))
      MOI.write_to_file(lp_file, "hard.lp")
    end

    if has_values(model)
      #if objective_value(model) <= cutoff
      if math_elim_mode in [2,3] && value.(model[:rank]) <= cutoff
        return true
      elseif math_elim_mode in [4,5] && value.(model[:W]) <= cutoff
        return true
      else
        return false
      end
    else # no values, return false (this is an assumption)
      return false
    end
  elseif status in [MOI.INFEASIBLE, MOI.INFEASIBLE_OR_UNBOUNDED]
    return false
  else
    error("The model was not solved correctly. Exiting with status $status.")
  end
  return false
end # checkMIP

"""
updateUsingMIPSolution!: Given a solution to the MIP, update best set of outcomes

Updates
  * best_outcomes
  * best_h2h
  * best_num_wins
  * best_rank
"""
function updateUsingMIPSolution!(model, k, t, schedule, h2h, W, num_playoff_teams,
    best_outcomes, best_h2h, best_num_wins, best_rank, math_elim_mode)
  status = termination_status(model)
  good_status = status in [MOI.OPTIMAL, MOI.OBJECTIVE_LIMIT, MOI.TIME_LIMIT]
  if !good_status || !has_values(model)
    error("The model was not solved correctly. Exiting with status $status.")
  end

  num_games_total = size(schedule, 1)
  num_teams = size(best_rank, 1)

  w = model[:w]
  for i = 1:num_teams
    best_num_wins[k, i] = Int(round(value.(w[i])))
  end

  if math_elim_mode in [2,4] # xit
    for game_ind = t:num_games_total
      i = schedule[game_ind,1]
      j = schedule[game_ind,2]
      xit = variable_by_name(model, "x[$i,$game_ind]")
      xjt = variable_by_name(model, "x[$j,$game_ind]")
      winner = isVal(value.(xit), 1) ? i : j

      best_outcomes[k, game_ind] = winner
    end
  elseif math_elim_mode in [3,5] # xij
    for i = 1:num_teams
      for j = 1:num_teams 
        if i == j
          continue
        end
        xij = variable_by_name(model, "x[$i,$j]")
        xji = variable_by_name(model, "x[$j,$i]")
        best_h2h[k,i,j] = Int(round(value.(xij)))
        best_h2h[k,j,i] = Int(round(value.(xji)))
      end
    end
  end

  if math_elim_mode in [2,3]
    best_rank[k] = Int(round(value.(model[:rank])))
  elseif math_elim_mode in [4,5]
    if (value.(model[:W]) <= W)
      best_rank[k] = num_playoff_teams
    else
      best_rank[k] = num_playoff_teams+1
    end
  end

end # updateUsingMIPSolution

"""
updateOtherUsingBestSolution!: Check whether other teams best schedule can be updated

Updates
  * best_outcomes
  * best_h2h
  * best_num_wins
  * best_rank
"""
function updateOthersUsingBestSolution!(k, t, schedule, num_playoff_teams,
    best_outcomes, best_h2h, best_num_wins, best_rank)
  num_games_total = length(schedule)
  num_teams = length(best_rank)

  num_wins = best_num_wins[k,:]
  W =  num_wins[num_playoff_teams]
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
      best_outcomes[i,:] = best_outcomes[k,:]
      best_h2h[i,:,:] = best_h2h[k,:,:]
      best_num_wins[i,:] = best_num_wins[k,:]
      best_rank[i] = rank_of_team[i]
    end
  end
end # updateOthersUsingBestSolution
