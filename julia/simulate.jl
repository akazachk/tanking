##################### 
# Tanking simulator #
#####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
include("utility.jl")
include("mathelim.jl")

DEBUG = false

"""
simulate: Simulates a season
  * num_teams
  * num_playoff_teams
  * num_rounds
  * num_replications
  * num_steps
  * gamma
  * breakpoint_list
  * nba_odds_list: list of odds to use by the NBA
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
"""
function simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, true_strength, mode, math_elim_mode=2)

  ## Set constants
  step_size                       = 1 / num_steps
  array_of_tanking_probabilities  = 0:step_size:1
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

  kend_out            = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # Kendall tau distance for bilevel ranking
  kend_nba_out        = zeros(Float64, num_steps+1, length(nba_odds_list), num_stats) # Kendall tau distance for NBA draft lottery system
  kend_gold_out       = zeros(Float64, 1, num_stats) # Kendall tau distance for Gold proposal
  kend_lenten_out     = zeros(Float64, 1, num_stats) # Kendall tau distance for Lenten proposal
  num_mips_out        = zeros(Float64, num_steps+1, num_stats) # number of MIPs solved in the course of a season
  games_tanked_out    = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # number of games tanked by each breakpoint
  already_tank_out    = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats) # number of teams tanking by each breakpoint
  math_eliminated_out = zeros(Float64, num_steps+1, num_games_total, num_stats) # number of teams mathematically eliminated by each game
  eff_eliminated_out  = zeros(Float64, num_steps+1, num_games_total, num_stats) # number of teams effectively eliminated by each game
  num_unelim_out      = zeros(Float64, num_steps+1, num_stats) # number teams eff elim then not not elim

  # Set min default (avg / stddev / max set to 0 is okay)
  BIG_NUMBER = max(num_teams^2, num_games_total)
  kend_out[:,:,min_stat]            = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
  kend_nba_out[:,:,min_stat]        = BIG_NUMBER * ones(num_steps+1, length(nba_odds_list))
  kend_gold_out[min_stat]           = BIG_NUMBER
  kend_lenten_out[min_stat]         = BIG_NUMBER
  games_tanked_out[:,:,min_stat]    = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
  already_tank_out[:,:,min_stat]    = BIG_NUMBER * ones(num_steps+1, length(breakpoint_list))
  math_eliminated_out[:,:,min_stat] = BIG_NUMBER * ones(num_steps+1, num_games_total)
  eff_eliminated_out[:,:,min_stat]  = BIG_NUMBER * ones(num_steps+1, num_games_total)
  num_mips_out[:,min_stat]          = BIG_NUMBER * ones(num_steps+1)
  num_unelim_out[:,min_stat]        = BIG_NUMBER * ones(num_steps+1)

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
  for step_ind in 1:length(array_of_tanking_probabilities)
    tank_perc = array_of_tanking_probabilities[step_ind]
    print("Simulating season with $tank_perc ratio of teams tanking\n")
    for rep = 1:num_replications
      print("\tReplication $rep/$num_replications, ratio $tank_perc\n")
      
      ## Set up stats for current replication
      for i = 1:num_teams
        stats[i,name_ind]       = i # team name
        stats[i,wins_ind]       = 0 # num wins
        stats[i,losses_ind]     = 0 # num losses 
        stats[i,games_left_ind] = num_team_games # num games left
        stats[i,elim_ind]       = -1 # when team is eliminated (game_ind)
        stats[i,win_pct_ind]    = 0.0 # win percentage
        if rand() > tank_perc
          stats[i,will_tank_ind] = 0 # will not tank (team is "moral")
        else
          stats[i,will_tank_ind] = 1 # will tank (team is "strategic")
        end
      end # set up stats array

      ## Prepare output for current replication
      num_eliminated      = 0
      num_math_elim       = 0
      num_eff_elim        = 0
      num_teams_tanking   = 0
      num_games_tanked    = 0
      outcome             = zeros(Int, num_games_total)
      is_eff_elim         = zeros(Bool, num_teams)
      is_math_elim        = zeros(Bool, num_teams)
      num_wins_since_elim = zeros(Int, num_teams) # for Gold ranking
      eff_elim_game       = -ones(Int, num_teams)
      math_elim_game      = -ones(Int, num_teams) # for Lenten ranking
      h2h                 = zeros(Int, num_teams, num_teams)
      h2h_left            = num_rounds * ones(Int, num_teams, num_teams)

      is_unelim           = zeros(Bool, num_teams)

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
      model = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode)
      # START DEBUG
      if DEBUG
        model2 = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, 2)
        model3 = setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, 3)
      end
      # END DEBUG

      ## Run one season
      for game_ind = 1:num_games_total
        # Current teams playing
        i = schedule[game_ind, 1] #games[round_ind, round_game_ind, 1]
        j = schedule[game_ind, 2] #games[round_ind, round_game_ind, 2]

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
          rank_of_team, team_in_pos = updateRank(stats, rank_of_team, team_in_pos, k, team_k_wins, num_teams, win_pct_ind, games_left_ind, h2h) # update rank
          
          # If team k wins and has been eliminated, increment num wins since elim (for Gold ranking)
          is_elim = math_elim_mode == 0 ? is_eff_elim[k] : is_math_elim[k] # note that we use math elimination if it has been calculated
          if team_k_wins && is_elim
            num_wins_since_elim[k] += 1
          end
        end
        #print("($i,$j) Team $i wins? $team_i_wins\n")
        #display([1:num_teams team_in_pos rank_of_team stats[team_in_pos,win_pct_ind]])

        # When the breakpoint for choosing a playoff ranking has been reached, set the ranking
        for r = 1:length(breakpoint_game_for_draft) 
          if game_ind == breakpoint_game_for_draft[r]
            draft_rank_of_team[:,r] = rank_of_team
            @views updateStats!(already_tank_out[step_ind, r, :], num_teams_tanking, num_replications)
            @views updateStats!(games_tanked_out[step_ind, r, :], num_games_tanked, num_replications)
          end
        end
        
        # Check whether teams are eliminated
        # Prepare for elimination calculations
        if math_elim_mode != 0
          # Update mathematical elimination
          updateHeuristicBestRank!(outcome[game_ind], game_ind, schedule, h2h, best_outcomes, best_h2h, best_num_wins, best_rank)
          fixOutcome!(model, outcome[game_ind], game_ind, schedule, math_elim_mode)
            ### START DEBUG
            if DEBUG
              fixOutcome!(model2, outcome[game_ind], game_ind, schedule, 2)
              fixOutcome!(model3, outcome[game_ind], game_ind, schedule, 3)
            end
            ### END DEBUG
        end
        
        # Find cutoff win percentage and set elimination status for teams eliminated after this game
        last_team = team_in_pos[num_playoff_teams]
        cutoff_avg = stats[last_team,win_pct_ind]

        # Loop through teams, checking each one's elimination status
        for k in 1:num_teams
          if math_elim_mode != 0 && !is_math_elim[k] && (math_elim_mode <= 3 || math_elim_mode >= 6 || k == 1)
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
                math_elim_mode, wins_ind, games_left_ind)
            num_mips += mips_used
            num_math_elim += is_math_elim[k]
            math_elim_game[k] = num_team_games - stats[k,games_left_ind] # number of games played at elimination
          end
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
                math_elim_mode, wins_ind, games_left_ind)
            num_mips += mips_used
            num_math_elim += is_math_elim[k]
        end
        
        #println("(step $step_ind, game $game_ind): num_mips: $num_mips\tnum_math_elim: $num_math_elim")
        
        # Update elimination stats
        @views updateStats!(math_eliminated_out[step_ind, game_ind, :], num_math_elim, num_replications)
        @views updateStats!(eff_eliminated_out[step_ind, game_ind, :], num_eff_elim, num_replications)
      end # iterate over games
      @views updateStats!(num_mips_out[step_ind,:], num_mips, num_replications)
      ## end of a season
      
      ## Make sure teams that were eliminated actually did not make the playoffs
      for i = 1:num_teams
        if (best_rank[i] < 0)
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
        curr_kend = kendtau_sorted(sorted_ranking[:,1], true_strength, mode)
        @views updateStats!(kend_out[step_ind, r, :], curr_kend, num_replications)

        # Now repeat with h2h
        #tmp_stats = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        #tmp_stats[:,1] = nonplayoff_teams_h2h
        #tmp_stats[:,2] = draft_rank_of_team_h2h[nonplayoff_teams, r]
        #sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
        #kend_out_h2h[step_ind, r] += kendtau_sorted(sorted_ranking[:,1], true_strength, mode) / num_replications
      end

      ## Compute the Kendall tau distance when the NBA randomization is used
      for r = 1:length(nba_odds_list)
        odds = copy(nba_odds_list[r])
        draft_order = runDraftLottery(nonplayoff_teams, odds, length(nonplayoff_teams))
        ranking_nba = draft_order[length(draft_order):-1:1]
        curr_kend = kendtau_sorted(ranking_nba, true_strength, mode)
        @views updateStats!(kend_nba_out[step_ind, r, :], curr_kend, num_replications)
      end

      ## When there is no tanking, compute the Gold and Lenten methods
      if (tank_perc == 0.0)
        ## Gold's ranking: based on number of wins since a team was eliminated
        ranking_gold = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        ranking_gold[:,1] = nonplayoff_teams
        ranking_gold[:,2] = -1 * num_wins_since_elim[nonplayoff_teams] # negative because teams with more wins need to be ranked worse (as they are given a _better_ draft pick)
        curr_kend = kendtau(ranking_gold, 2, true_strength, mode)
        @views updateStats!(kend_gold_out, curr_kend, num_replications)

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
        @views updateStats!(kend_lenten_out, curr_kend, num_replications)
      end # check if tank_perc == 0 (for computing Lenten and Gold rankings)
    end # do replications

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

    num_mips_out[step_ind, stddev_stat]   -= num_mips_out[step_ind, avg_stat]^2
    num_unelim_out[step_ind, stddev_stat] -= num_unelim_out[step_ind, avg_stat]^2

    if (tank_perc == 0.0)
      kend_gold_out[stddev_stat]   -= kend_gold_out[avg_stat]^2
      kend_lenten_out[stddev_stat] -= kend_lenten_out[avg_stat]^2
    end
    #println("num_mips_out: ", num_mips_out[step_ind,:])
    #quit()
  end # looping over tanking percentages

  return kend_out, kend_nba_out, kend_gold_out, kend_lenten_out,
      games_tanked_out, already_tank_out, 
      math_eliminated_out, eff_eliminated_out, 
      num_mips_out, num_unelim_out
end # simulate
