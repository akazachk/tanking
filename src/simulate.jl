##################### 
# Tanking simulator #
#####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
#include("utility.jl")
#include("mathelim.jl")

import Random.randperm

#import Gurobi
#const GRB_ENV = Gurobi.Env()

DEBUG = false

"""
simulate: Simulates a season

Teams play each other in rounds, consisting of each team playing every other team
Each team may or may not tank at all (decided by a tanking percentage)

If a team might tank, it will tank when it decides it has no chance of making the playoffs
(after a minimum number of games have been played, currently set to half their games)
This is when it would not be enough for the team to win all its remaining games 
to have a win percentage at least as good as the last playoff team
(assuming that cutoff win percentage remains the same)

Assumptions:
1. No simultaneous games
2. Teams keep same relative true ranking throughout season
3. No conferences / divisions
4. No home / away games (i.e., no home / away advantages)

Parameters
---
  * num_teams
  * num_playoff_teams
  * num_rounds
  * num_replications
  * num_steps
  * gamma
  * breakpoint_list
  * nba_odds_list: list of odds to use by the NBA
  * nba_num_lottery: number of teams for which the lottery decides draft position
  * true_strength
  * mode: defines how the ranking and winner determination works
    1 or 2: 
      When two non-tanking teams or two tanking teams play each other, 
      the better team wins with probability gamma
      When a tanking team plays a non-tanking team, the tanking team always loses
    3 or 4:
      Variants of (Zermelo-)Bradley-Terry model used to determine who wins each game
  * math_elim_mode: 
    0: use effective elimination, 
    1: use mathematical elimination, but calculated by heuristics only, 
    2: use math elim, binary MIP, team-wise formulation
    3: use math elim, general integer MIP, team-wise formulation
    4: use math elim, binary MIP, cutoff formulation
    5: use math elim, general integer MIP, cutoff formulation
    <0: use effective elimination, but calculate mathematical elimination
  * selected_steps: which steps to do
  * env: Gurobi environment
  * ONLY_RETURN_WIN_PCT: do not calculate all of the return values, only avg_win_pct

Returns
---
  * kend
  * kend_nba
  * kend_gold
  * kend_lenten
  * games_tanked
  * already_tank
  * math_eliminated
  * eff_eliminated
  * num_mips
  * num_unelim
  * avg_rank_strat
  * avg_rank_moral
  * avg_elim_rank_strat
  * avg_elim_rank_moral
  * avg_diff_rank_strat
  * avg_diff_rank_moral
  * num_missing_case
"""
function simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, nba_num_lottery, true_strength, mode, math_elim_mode=-2, selected_steps=nothing, env=nothing, ONLY_RETURN_WIN_PCT=false)

  ## Set constants
  step_size                       = 1 / num_steps
  array_of_tanking_probabilities  = 0:step_size:1
  if isa(selected_steps, Number)
    selected_steps = [selected_steps]
  elseif !isa(selected_steps,Array) && !isa(selected_steps,UnitRange)
    selected_steps = 1:length(array_of_tanking_probabilities)
  end
  num_games_per_round             = Int(num_teams * (num_teams - 1) / 2)
  num_games_total                 = num_rounds * num_games_per_round
  num_team_games                  = num_rounds * (num_teams - 1)
  max_games_remaining             = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen

  breakpoint_game_for_draft = [round(breakpoint_list[i] * num_games_total) for i in 1:length(breakpoint_list)]

  ## Prepare output
  # For each stat, keep: 1. avg, 2. stddev, 3. min, 4. max
  num_stats   = 4
  avg_stat    = 1
  stddev_stat = 2
  min_stat    = 3
  max_stat    = 4
  BIG_NUMBER = max(num_teams^2, num_games_total)

  if ONLY_RETURN_WIN_PCT
    win_pct_out = zeros(Float64, num_steps+1, num_teams, num_stats)
    win_pct_out[:,:,min_stat] = BIG_NUMBER * ones(num_steps+1, num_teams)
  else
    kend_out            = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # Kendall tau distance for bilevel ranking
    kend_nba_out        = zeros(Float64, num_steps+1, length(nba_odds_list), num_stats) # Kendall tau distance for NBA draft lottery system
    kend_gold_out       = zeros(Float64, 1, num_stats) # Kendall tau distance for Gold proposal
    kend_lenten_out     = zeros(Float64, num_steps+1, num_stats) # Kendall tau distance for Lenten proposal
    num_mips_out        = zeros(Float64, num_steps+1, num_stats) # number of MIPs solved in the course of a season
    games_tanked_out    = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # number of games tanked by each breakpoint
    already_tank_out    = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # number of teams tanking by each breakpoint
    math_eliminated_out = zeros(Float64, num_steps+1, num_games_total, num_stats) # number of teams mathematically eliminated by each game
    eff_eliminated_out  = zeros(Float64, num_steps+1, num_games_total, num_stats) # number of teams effectively eliminated by each game
    num_unelim_out      = zeros(Float64, num_steps+1, num_stats) # number teams eff elim then not not elim
    avg_rank_strat_out  = zeros(Float64, num_steps+1, num_stats) # avg rank of strategic teams
    avg_rank_moral_out  = zeros(Float64, num_steps+1, num_stats) # avg rank of moral teams
    avg_elim_rank_strat_out = zeros(Float64, num_steps+1, num_stats) # avg rank of eliminated strategic teams
    avg_elim_rank_moral_out = zeros(Float64, num_steps+1, num_stats) # avg rank of eliminated moral teams
    avg_diff_rank_strat_out = zeros(Float64, num_steps+1, num_stats) # avg of diff between true and calc ranks of strategic teams
    avg_diff_rank_moral_out = zeros(Float64, num_steps+1, num_stats) # avg of diff between true and calc ranks of moral teams
    num_missing_case_out = zeros(Float64, num_steps+1, num_stats) # number teams eff elim then not not elim
    SHOULD_KEEP_H2H_OUT = false
    if SHOULD_KEEP_H2H_OUT
      h2h_out              = zeros(Float64, num_steps+1, num_teams, num_teams, num_stats) # average head-to-head statistics
    end

    # Set min default (avg / stddev / max set to 0 is okay)
    kend_out[:,:,min_stat]            = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
    kend_nba_out[:,:,min_stat]        = BIG_NUMBER * ones(num_steps+1, length(nba_odds_list))
    kend_gold_out[min_stat]           = BIG_NUMBER
    kend_lenten_out[:,min_stat]       = BIG_NUMBER * ones(num_steps+1)
    games_tanked_out[:,:,min_stat]    = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
    already_tank_out[:,:,min_stat]    = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
    math_eliminated_out[:,:,min_stat] = BIG_NUMBER * ones(num_steps+1, num_games_total)
    eff_eliminated_out[:,:,min_stat]  = BIG_NUMBER * ones(num_steps+1, num_games_total)
    num_mips_out[:,min_stat]          = BIG_NUMBER * ones(num_steps+1)
    num_unelim_out[:,min_stat]        = BIG_NUMBER * ones(num_steps+1)
    avg_rank_strat_out[:,min_stat]    = BIG_NUMBER * ones(num_steps+1)
    avg_rank_moral_out[:,min_stat]    = BIG_NUMBER * ones(num_steps+1)
    avg_elim_rank_strat_out[:,min_stat] = BIG_NUMBER * ones(num_steps+1)
    avg_elim_rank_moral_out[:,min_stat] = BIG_NUMBER * ones(num_steps+1)
    avg_diff_rank_strat_out[:,min_stat] = BIG_NUMBER * ones(num_steps+1)
    avg_diff_rank_moral_out[:,min_stat] = BIG_NUMBER * ones(num_steps+1)
    if SHOULD_KEEP_H2H_OUT
      h2h_out[:,:,:,min_stat]           = BIG_NUMBER * ones(num_steps+1, num_teams, num_teams)
    end
  end # ONLY_RETURN_WIN_PCT

  ## Set up for game order 
  # Alternative to below is using Combinatorics; ord_games = collect(combinations(1:num_teams,2))
  schedule = Matrix{Int64}(undef, num_rounds * num_games_per_round, 2) # schedule of games
  ord_games = Matrix{Int64}(undef, num_games_per_round, 2) # all the games played per round [internal]
  curr_game_ind = 1 
  for i in 1:num_teams
    for j in (i+1):num_teams
      ord_games[curr_game_ind,1] = i
      ord_games[curr_game_ind,2] = j
      curr_game_ind = curr_game_ind + 1
    end
  end

  ## For ranking teams
  # stats[:,1] is team "name" (index of the team)
  # stats[:,2] is num wins
  # stats[:,3] is num losses
  # stats[:,4] is num games remaining
  # stats[:,5] is when team is eliminated 
  # stats[:,6] is win percentage
  # stats[:,7] is indicator for whether team tanks
  size_of_stats   = 7
  name_ind        = 1
  wins_ind        = 2
  losses_ind      = 3
  games_left_ind  = 4
  elim_ind        = 5
  win_pct_ind     = 6
  will_tank_ind   = 7
  stats = Matrix{Any}(undef, num_teams, size_of_stats)
  draft_rank_of_team = Array{Any}(undef, num_teams, length(breakpoint_game_for_draft))

  ## Begin calculations
  decide_tanking_with_prob = true
  if num_steps == num_teams
    decide_tanking_with_prob = false
  end
  for step_ind in selected_steps
    tank_perc = array_of_tanking_probabilities[step_ind]
    num_repl_for_avg = 0
    print("Simulating season with $tank_perc ratio of teams tanking\n")
    for rep = 1:num_replications
      #print("\tReplication $rep/$num_replications, ratio $tank_perc\n")
      
      ## Prepare to decide strategic teams
      num_tanking = 0
      teamperm = 0
      if !decide_tanking_with_prob
        num_tanking = Int(round(tank_perc * num_teams))
        teamperm = randperm(num_teams)
        teamperm = teamperm[1:num_tanking]
      end

      ## Set up stats for current replication
      for i = 1:num_teams
        stats[i,name_ind]       = i # team name
        stats[i,wins_ind]       = 0 # num wins
        stats[i,losses_ind]     = 0 # num losses 
        stats[i,games_left_ind] = num_team_games # num games left
        stats[i,elim_ind]       = -1 # when team is eliminated (game_ind)
        stats[i,win_pct_ind]    = 0.0 # win percentage
        if decide_tanking_with_prob
          if rand() > tank_perc
            stats[i,will_tank_ind] = 0 # will not tank (team is "moral")
          else
            stats[i,will_tank_ind] = 1 # will tank (team is "strategic")
            num_tanking += 1
          end
        else
          if num_tanking > 0 && i in teamperm
            stats[i,will_tank_ind] = 1 # will tank (team is "strategic")
          else
            stats[i,will_tank_ind] = 0 # will not tank (team is "moral")
          end
        end
      end # set up stats array

      ## Prepare output for current replication
      num_eliminated      = 0
      num_math_elim       = 0
      num_eff_elim        = 0
      num_teams_tanking   = 0
      num_games_tanked    = 0
      rank_strat          = 0
      rank_moral          = 0
      elim_rank_strat     = 0
      elim_rank_moral     = 0
      diff_rank_strat     = 0
      diff_rank_moral     = 0
      outcome             = zeros(Int, num_games_total)
      is_eff_elim         = zeros(Bool, num_teams)
      is_math_elim        = zeros(Bool, num_teams)
      num_wins_since_elim = zeros(Int, num_teams) # for Gold ranking
      eff_elim_game       = -ones(Int, num_teams)
      math_elim_game      = -ones(Int, num_teams) # for Lenten ranking
      h2h                 = zeros(Int, num_teams, num_teams)
      h2h_left            = num_rounds * ones(Int, num_teams, num_teams)

      is_unelim           = zeros(Bool, num_teams) # when effective elimination is wrong
      math_elim_game_total = -ones(Int, num_teams) # game (within whole schedule) that team i is eliminated
      num_missing_case    = 0

      ## Set up initial ranking
      rank_of_team = sortperm(randn(num_teams)) # initial ranking (returns rank of team i)
      team_in_pos = Array{Int}(undef, num_teams) # inverse ranking (returns team that is in position i)
      for i = 1:num_teams
        team_in_pos[rank_of_team[i]] = i
      end

      ## Set random game order for this replication
      for round_ind = 1:num_rounds
        perm = sortperm(randn(num_games_per_round))
        start_round = num_games_per_round * (round_ind - 1) + 1
        end_round = num_games_per_round * round_ind
        schedule[start_round:end_round, :] = ord_games[perm,:]
      end

      ## Set up MIP
      best_outcomes = zeros(Int, num_teams, num_games_total)
      best_h2h      = zeros(Int, num_teams, num_teams, num_teams)
      best_num_wins = zeros(Int, num_teams, num_teams)
      best_rank     = zeros(Int, num_teams)
      num_mips = 0
      model = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, abs(math_elim_mode), env)
      # START DEBUG
      if DEBUG
        model2 = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, 2, env)
        model3 = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, 3, env)
      end
      # END DEBUG

      ## Run one season
      for game_ind = 1:num_games_total
        # Current teams playing
        i = schedule[game_ind, 1] #games[round_ind, round_game_ind, 1]
        j = schedule[game_ind, 2] #games[round_ind, round_game_ind, 2]

        # Decide if incomplete case has been encountered:
        # (1) it is after game delta
        # (2) team i has been eliminated
        # (3) team i was not eliminated at the breakpoint game delta
        # (4) team j was ranked worse than team i at the breakpoint game delta
        # (5) team j is a contender (has neither been eliminated nor is guaranteed to make the playoffs)
        # (6) some remaining contender was ranked better than team i at the breakpoint game delta
        if math_elim_mode != 0
          for r = 1:length(breakpoint_game_for_draft)
            delta = breakpoint_game_for_draft[r]
            # Check condition (1)
            if game_ind <= delta
              continue
            end

            # Check conditions (2)-(5) for team i
            cond23 = (math_elim_game_total[i] > delta) # if math_elim_game_total[i] < 0, i.e., the team is not elim, this is false
            cond4 = (draft_rank_of_team[i,r] < draft_rank_of_team[j,r])
            cond5 = teamIsContender(j, game_ind, schedule, stats, outcome, h2h, math_elim_game_total, num_playoff_teams, wins_ind, games_left_ind)
            if cond23 && cond4 && cond5
              # Check condition (6) that some remaining contender was ranked better than team i at the breakpoint game delta
              exists_contender = false
              for k = 1:num_teams
                if (draft_rank_of_team[i,r] > draft_rank_of_team[k,r]) && teamIsContender(k, game_ind, schedule, stats, outcome, h2h, math_elim_game_total, num_playoff_teams, wins_ind, games_left_ind) 
                  exists_contender = true
                  break
                end
              end # loop over teams looking for contender for condition (6)

              num_missing_case += exists_contender
            end # check if conditions (2)-(5) are met for team i

            # Check condtions (2)-(5) for team j
            cond23 = (math_elim_game_total[j] > delta) # if math_elim_game_total[i] < 0, i.e., the team is not elim, this is false
            cond4 = (draft_rank_of_team[j,r] < draft_rank_of_team[i,r])
            cond5 = teamIsContender(i, game_ind, schedule, stats, outcome, h2h, math_elim_game_total, num_playoff_teams, wins_ind, games_left_ind)
            if cond23 && cond4 && cond5
              # Check condition (6) that some remaining contender was ranked better than team j at the breakpoint game delta
              exists_contender = false
              for k = 1:num_teams
                if (draft_rank_of_team[j,r] > draft_rank_of_team[k,r]) && teamIsContender(k, game_ind, schedule, stats, outcome, h2h, math_elim_game_total, num_playoff_teams, wins_ind, games_left_ind) 
                  exists_contender = true
                  break
                end
              end # loop over teams looking for contender for condition (6)

              num_missing_case += exists_contender
            end # check if conditions (2)-(5) are met for team j
          end # loop over breakpoint games
        end # ensure math_elim_mode != 0

        # Decide who wins the game
        team_i_wins = teamWillWin(i, j, stats, gamma, true_strength, mode, elim_ind, will_tank_ind)

        # Check if any teams tanked this game
        team_i_is_tanking = teamIsTanking(i, stats, elim_ind, will_tank_ind)
        team_j_is_tanking = teamIsTanking(j, stats, elim_ind, will_tank_ind)
        if team_i_is_tanking || team_j_is_tanking
          num_games_tanked += 1
        end

        # Do updates
        outcome[game_ind] = team_i_wins ? i : j
        h2h[i,j]         += team_i_wins
        h2h[j,i]         += !team_i_wins
        h2h_left[i,j]    -= 1
        h2h_left[j,i]    -= 1
        for k in [i,j]
          team_k_wins               = (k == i) ? team_i_wins : !team_i_wins
          stats[k,wins_ind]        += team_k_wins
          stats[k,losses_ind]      += !team_k_wins
          stats[k,games_left_ind]  -= 1 # one fewer game remaining
          stats[k,win_pct_ind]      = stats[k,wins_ind] / (num_team_games - stats[k,games_left_ind]) # update current win pct
          rank_of_team, team_in_pos = updateRank(rank_of_team, team_in_pos, stats, k, team_k_wins, num_teams, win_pct_ind, games_left_ind, h2h) # update rank
          
          # If team k wins and has been eliminated, increment num wins since elim (for Gold ranking)
          # NB: this is NOT the same as the is_elim used for tanking decisions later
          # Namely, we use math elimination here if it has been calculated, even if math_elim_mode < 0
          is_elim = math_elim_mode == 0 ? is_eff_elim[k] : is_math_elim[k]
          if team_k_wins && is_elim
            num_wins_since_elim[k] += 1
          end
        end
        #print("($i,$j) Team $i wins? $team_i_wins\n")
        #display([1:num_teams team_in_pos rank_of_team stats[team_in_pos,win_pct_ind]])

        # Check whether teams are eliminated
        # Prepare for elimination calculations
        if math_elim_mode != 0
          # Update mathematical elimination
          updateHeuristicBestRank!(outcome[game_ind], game_ind, schedule, h2h, best_outcomes, best_h2h, best_num_wins, best_rank)
          fixOutcome!(model, outcome[game_ind], game_ind, schedule, abs(math_elim_mode))
            ### START DEBUG
            if DEBUG
              fixOutcome!(model2, outcome[game_ind], game_ind, schedule, 2)
              fixOutcome!(model3, outcome[game_ind], game_ind, schedule, 3)
            end
            ### END DEBUG
        end
        
        # Find cutoff win percentage and set elimination status for teams eliminated after this game
        last_team = num_playoff_teams > 0 ? team_in_pos[num_playoff_teams] : -1
        cutoff_avg = num_playoff_teams > 0 ? stats[last_team,win_pct_ind] : 1.1

        # Loop through teams, checking each one's elimination status
        for k in 1:num_teams
          if math_elim_mode != 0 && !is_math_elim[k] && (abs(math_elim_mode) <= 3 || abs(math_elim_mode) >= 6 || k == 1)
            ### START DEBUG
            if DEBUG
              best_outcomes2 = copy(best_outcomes)
              best_h2h2 = copy(best_h2h)
              best_num_wins2 = copy(best_num_wins)
              best_rank2 = copy(best_rank)
              (is_elim2, mips_used2) = teamIsMathematicallyEliminated!(k, game_ind, 
                  schedule, stats, outcome, h2h, best_outcomes2, best_h2h2, best_num_wins2, best_rank2, model2,
                  num_teams, num_playoff_teams, num_team_games, num_games_total,
                  2, wins_ind, games_left_ind)

              best_outcomes3 = copy(best_outcomes)
              best_h2h3 = copy(best_h2h)
              best_num_wins3 = copy(best_num_wins)
              best_rank3 = copy(best_rank)
              (is_elim3, mips_used3) = teamIsMathematicallyEliminated!(k, game_ind, 
                  schedule, stats, outcome, h2h, best_outcomes3, best_h2h3, best_num_wins3, best_rank3, model3,
                  num_teams, num_playoff_teams, num_team_games, num_games_total,
                  3, wins_ind, games_left_ind)

              if (is_elim2 != is_elim3)
                println("step $step_ind, game $game_ind, team $k: is_elim2: $is_elim2, is_elim3: $is_elim3")
                ## Team k is eliminated in one but not the other model; impossible. We should be able to find the schedule in which the team is not eliminated
                println("using 2: best_num_wins = ", best_num_wins2[k, k], ", best_rank2 = ", best_rank2[k])
                println("using 3: best_num_wins = ", best_num_wins3[k, k], ", best_rank3 = ", best_rank3[k])

                if !is_elim2
                  # Check team k's rank
                end
              end
              @assert(is_elim2 == is_elim3)
            end
            ### END DEBUG

            # Mathematical elimination: solve MIP if heuristic does not find good schedule
            (is_math_elim[k], mips_used) = teamIsMathematicallyEliminated!(k, game_ind, 
                schedule, stats, outcome, h2h, best_outcomes, best_h2h, best_num_wins, best_rank, model,
                num_teams, num_playoff_teams, num_team_games, num_games_total,
                abs(math_elim_mode), wins_ind, games_left_ind)
            num_mips += mips_used
            if is_math_elim[k]
              num_math_elim += 1
              math_elim_game[k] = num_team_games - stats[k,games_left_ind] # number of games played at elimination
              math_elim_game_total[k] = game_ind
            end
          end # check that team is not yet math elim (and that we should check math elim)
          if !is_eff_elim[k]
            is_eff_elim[k] = teamIsEffectivelyEliminated(stats[k,wins_ind], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
            num_eff_elim += is_eff_elim[k]
            eff_elim_game[k] = num_team_games - stats[k,games_left_ind] # number of games played at elimination
          end

          # If current team is eliminated, and it has not been recorded before, do so
          is_elim = math_elim_mode > 0 ? is_math_elim[k] : is_eff_elim[k]
          if is_elim && stats[k,elim_ind] < 0
            num_eliminated += 1
            stats[k,elim_ind] = num_team_games - stats[k,games_left_ind] # number of games played at elimination
            num_teams_tanking += stats[k,will_tank_ind] == 1
          end

          # Check whether an effectively eliminated team is no longer eliminated
          if is_eff_elim[k] && !is_unelim[k] && (rank_of_team[k] <= num_playoff_teams)
            is_unelim[k] = true
          end
        end # check elimination

        # For math elimination modes that use the cutoff formulation, we solve one MIP per game
        if false
            (is_math_elim[k], mips_used) = teamIsMathematicallyEliminated!(k, game_ind, 
                schedule, stats, outcome, h2h, best_outcomes, best_h2h, best_num_wins, best_rank, model,
                num_teams, num_playoff_teams, num_team_games, num_games_total,
                abs(math_elim_mode), wins_ind, games_left_ind)
            num_mips += mips_used
            num_math_elim += is_math_elim[k]
        end
        
        #println("(step $step_ind, game $game_ind): num_mips: $num_mips\tnum_math_elim: $num_math_elim")

        if ONLY_RETURN_WIN_PCT
          continue
        end

        # When the breakpoint for choosing a playoff ranking has been reached, set the ranking
        for r = 1:length(breakpoint_game_for_draft) 
          if game_ind == breakpoint_game_for_draft[r]
            draft_rank_of_team[:,r] = rank_of_team
            @views updateStats!(already_tank_out[step_ind, r, :], num_teams_tanking, num_replications)
            @views updateStats!(games_tanked_out[step_ind, r, :], num_games_tanked, num_replications)
          end
        end
        
        # Update elimination stats
        @views updateStats!(math_eliminated_out[step_ind, game_ind, :], num_math_elim, num_replications)
        @views updateStats!(eff_eliminated_out[step_ind, game_ind, :], num_eff_elim, num_replications)
      end # iterate over games
      ## end of a season

      if ONLY_RETURN_WIN_PCT
        for i in 1:num_teams
          # Make sure we use sorted ranking
          @views updateStats!(win_pct_out[step_ind, i, :], 100 * stats[team_in_pos[i], win_pct_ind], num_replications)
        end
        continue
      end

      ## Update num_mips, num_missing_case, and h2h stats
      @views updateStats!(num_mips_out[step_ind,:], num_mips, num_replications)
      @views updateStats!(num_missing_case_out[step_ind, :], num_missing_case, num_replications)
      if SHOULD_KEEP_H2H_OUT
        for i = 1:num_teams
          for j = i+1:num_teams
            @views updateStats!(h2h_out[step_ind, i, j, :], h2h[i,j], num_replications)
            @views updateStats!(h2h_out[step_ind, j, i, :], h2h[j,i], num_replications)
            @assert( h2h[i,j] + h2h[j,i] == num_rounds )
          end
        end
      end
      
      ## Make sure teams that were eliminated actually did not make the playoffs
      for i = 1:num_teams
        if (best_rank[i] < 0)
          if (rank_of_team[i] <= num_playoff_teams)
            println("violation found for team $i")
            println("best_rank: ", best_rank)
            println("rank_of_team: ", rank_of_team)
          end
          @assert(rank_of_team[i] > num_playoff_teams)
        end
      end

      ## Update uneliminated stats
      num_unelim = sum(is_unelim)
      @views updateStats!(num_unelim_out[step_ind,:], num_unelim, num_replications)
          
      ## Get non-playoff teams at end of season
      nonplayoff_teams = team_in_pos[num_playoff_teams+1:num_teams]
      #nonplayoff_teams_h2h = team_in_pos_h2h[num_playoff_teams+1:num_teams]
      for r = 1:length(breakpoint_game_for_draft)
        tmp_stats = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        tmp_stats[:,1] = nonplayoff_teams
        tmp_stats[:,2] = draft_rank_of_team[nonplayoff_teams, r]
        sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
        curr_kend = kendtau_sorted(sorted_ranking[:,1], true_strength, mode=mode)
        @views updateStats!(kend_out[step_ind, r, :], curr_kend, num_replications)
        ## DEBUG
        #println("wins = ", stats[:, wins_ind], "\tkend = $curr_kend\tavg_kend = ", kend_out[step_ind, r, avg_stat])
        ## DEBUG

        # Now repeat with h2h
        #tmp_stats = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        #tmp_stats[:,1] = nonplayoff_teams_h2h
        #tmp_stats[:,2] = draft_rank_of_team_h2h[nonplayoff_teams, r]
        #sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
        #kend_out_h2h[step_ind, r] += kendtau_sorted(sorted_ranking[:,1], true_strength, mode=mode) / num_replications
      end

      ## Compute the Kendall tau distance when the NBA randomization is used
      for r = 1:length(nba_odds_list)
        odds = copy(nba_odds_list[r])
        draft_order = runDraftLottery(nonplayoff_teams, odds, length(nonplayoff_teams), nba_num_lottery[r])
        ranking_nba = draft_order[length(draft_order):-1:1]
        curr_kend = kendtau_sorted(ranking_nba, true_strength, mode=mode)
        @views updateStats!(kend_nba_out[step_ind, r, :], curr_kend, num_replications)
      end

      ## Compute Lenten Kendall tau distance
      ## Lenten ranking: the fewer games a team has played when it is eliminated, the higher it is in the draft
      ranking_lenten = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
      ranking_lenten[:,1] = nonplayoff_teams
      ranking_lenten[:,2] = math_elim_mode == 0 ? eff_elim_game[nonplayoff_teams] : math_elim_game[nonplayoff_teams] # not using negative, because we will sort high-to-low later
      for tmp_team_ind = 1:num_teams - num_playoff_teams
        if ranking_lenten[tmp_team_ind,2] < 0 # was never eliminated
          # Team was not eliminated but did not make the playoffs
          # It is so far unranked from Lenten perspective
          # Should be ranked higher than teams eliminated earlier
          ranking_lenten[tmp_team_ind,2] = num_games_total + 1
        end
      end
      curr_kend = kendtau(ranking_lenten, 2, true_strength, mode) # having a higher elimination index means Lenten ranks the team better (since it was eliminated later), i.e., it has a worse draft pick
      @views updateStats!(kend_lenten_out[step_ind,:], curr_kend, num_replications)

      ## When there is no tanking, compute the Gold and Lenten methods
      if (tank_perc == 0.0)
        ## Gold's ranking: based on number of wins since a team was eliminated
        ranking_gold = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        ranking_gold[:,1] = nonplayoff_teams
        ranking_gold[:,2] = -1 * num_wins_since_elim[nonplayoff_teams] # negative because teams with more wins need to be ranked worse (as they are given a _better_ draft pick)
        curr_kend = kendtau(ranking_gold, 2, true_strength, mode)
        @views updateStats!(kend_gold_out, curr_kend, num_replications)
      end # check if tank_perc == 0 (for computing Lenten and Gold rankings)

      ## Compute true vs calculated ranking among nonplayoff teams
      # We look at eliminated teams, and we recalculated the number of eliminated teams that tank
      # This number may be different from num_teams_tanking because we want to account for cases
      # in which teams would have tanked if given the chance, but did not
      # (otherwise, we are overstating the effect that tanking teams can have)
      num_elim_tanking = 0
      for i in 1:num_teams
        if (stats[i,will_tank_ind] == 1)
          rank_strat += rank_of_team[i] / num_tanking
        else
          rank_moral += rank_of_team[i] / (num_teams - num_tanking)
        end
      end
      for i in nonplayoff_teams
        if (stats[i,will_tank_ind] == 1)
          num_elim_tanking += 1
          elim_rank_strat += rank_of_team[i]
          diff_rank_strat += (rank_of_team[i] - i)
        else
          elim_rank_moral += rank_of_team[i]
          diff_rank_moral += (rank_of_team[i] - i)
        end
      end
      if num_elim_tanking > 0
        elim_rank_strat /= num_elim_tanking
        diff_rank_strat /= num_elim_tanking
      end
      if num_teams - num_playoff_teams - num_elim_tanking > 0
        elim_rank_moral /= (num_teams - num_playoff_teams - num_elim_tanking)
        diff_rank_moral /= (num_teams - num_playoff_teams - num_elim_tanking)
      end

      @views updateStats!(avg_rank_strat_out[step_ind,:], rank_strat, num_replications)
      @views updateStats!(avg_rank_moral_out[step_ind,:], rank_moral, num_replications)
      if (num_elim_tanking > 0) && (num_teams - num_playoff_teams - num_elim_tanking > 0)
        num_repl_for_avg += 1
        @views updateStats!(avg_elim_rank_strat_out[step_ind,:], elim_rank_strat, 1)
        @views updateStats!(avg_elim_rank_moral_out[step_ind,:], elim_rank_moral, 1)
        @views updateStats!(avg_diff_rank_strat_out[step_ind,:], diff_rank_strat, 1)
        @views updateStats!(avg_diff_rank_moral_out[step_ind,:], diff_rank_moral, 1)
      end
    end # do replications

    if ONLY_RETURN_WIN_PCT
      for i in 1:num_teams
        win_pct_out[step_ind, i, stddev_stat] -= win_pct_out[step_ind, i, avg_stat]^2
      end
      continue
    end

    if num_repl_for_avg > 0
      ## Update average for those stats that were not over all replications
      # avg
      avg_elim_rank_strat_out[step_ind,avg_stat] /= num_repl_for_avg
      avg_elim_rank_moral_out[step_ind,avg_stat] /= num_repl_for_avg
      avg_diff_rank_strat_out[step_ind,avg_stat] /= num_repl_for_avg
      avg_diff_rank_moral_out[step_ind,avg_stat] /= num_repl_for_avg
      # stddev
      avg_elim_rank_strat_out[step_ind,stddev_stat] /= num_repl_for_avg
      avg_elim_rank_moral_out[step_ind,stddev_stat] /= num_repl_for_avg
      avg_diff_rank_strat_out[step_ind,stddev_stat] /= num_repl_for_avg
      avg_diff_rank_moral_out[step_ind,stddev_stat] /= num_repl_for_avg
    end

    ## Compute standard deviation
    # Var[X] = E[X^2] - E[X]^2 
    # = (1/n) * (sum x_i^2) - ((1/n) * (sum x_i))^2
    # StdDev[X] = sqrt(Var[X])
    for r = 1:length(breakpoint_game_for_draft)
      kend_out[step_ind, r, stddev_stat]            -= kend_out[step_ind, r, avg_stat]^2
      games_tanked_out[step_ind, r, stddev_stat]    -= games_tanked_out[step_ind, r, avg_stat]^2
      already_tank_out[step_ind, r, stddev_stat]    -= already_tank_out[step_ind, r, avg_stat]^2
      math_eliminated_out[step_ind, r, stddev_stat] -= math_eliminated_out[step_ind, r, avg_stat]^2
      eff_eliminated_out[step_ind, r, stddev_stat]  -= eff_eliminated_out[step_ind, r, avg_stat]^2
    end

    for r = 1:length(nba_odds_list)
      kend_nba_out[step_ind, r, stddev_stat] -= kend_nba_out[step_ind, r, avg_stat]^2
    end

    num_mips_out[step_ind, stddev_stat] -= num_mips_out[step_ind, avg_stat]^2
    num_unelim_out[step_ind, stddev_stat] -= num_unelim_out[step_ind, avg_stat]^2
    kend_lenten_out[step_ind, stddev_stat] -= kend_lenten_out[step_ind, avg_stat]^2
    avg_rank_strat_out[step_ind, stddev_stat] -= avg_rank_strat_out[step_ind, avg_stat]^2
    avg_rank_moral_out[step_ind, stddev_stat] -= avg_rank_moral_out[step_ind, avg_stat]^2
    avg_elim_rank_strat_out[step_ind, stddev_stat] -= avg_elim_rank_strat_out[step_ind, avg_stat]^2
    avg_elim_rank_moral_out[step_ind, stddev_stat] -= avg_elim_rank_moral_out[step_ind, avg_stat]^2
    avg_diff_rank_strat_out[step_ind, stddev_stat] -= avg_diff_rank_strat_out[step_ind, avg_stat]^2
    avg_diff_rank_moral_out[step_ind, stddev_stat] -= avg_diff_rank_moral_out[step_ind, avg_stat]^2
    num_missing_case_out[step_ind, stddev_stat] -= num_missing_case_out[step_ind, avg_stat]^2
    if SHOULD_KEEP_H2H_OUT
      for i = 1:num_teams
        for j = i+1:num_teams
          h2h_out[step_ind,i,j,stddev_stat] -= h2h_out[step_ind,i,j,avg_stat]^2
          h2h_out[step_ind,j,i,stddev_stat] -= h2h_out[step_ind,j,i,avg_stat]^2
        end
      end
    end

    if (tank_perc == 0.0)
      kend_gold_out[stddev_stat] -= kend_gold_out[avg_stat]^2
    end
  end # looping over tanking percentages
  ### DEBUG
  #if SHOULD_KEEP_H2H_OUT
  #  println("h2h_out = ", h2h_out[1, :, :, avg_stat])
  #end
  ### DEBUG

  if ONLY_RETURN_WIN_PCT
    return win_pct_out
  end

  return kend_out, kend_nba_out, kend_gold_out, kend_lenten_out,
      games_tanked_out, already_tank_out, 
      math_eliminated_out, eff_eliminated_out, 
      num_mips_out, num_unelim_out,
      avg_rank_strat_out, avg_rank_moral_out,
      avg_elim_rank_strat_out, avg_elim_rank_moral_out,
      avg_diff_rank_strat_out, avg_diff_rank_moral_out,
      num_missing_case_out
end # simulate
