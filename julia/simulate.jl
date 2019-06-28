##################### 
# Tanking simulator #
#####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
include("utility.jl")
include("mathelim.jl")

function simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, true_strength, mode, return_h2h=false)
  ###
  # Simulates a season
  #
  # Teams play each other in rounds, consisting of each team playing every other team
  # Each team may or may not tank at all (decided by a tanking percentage)
  #
  # If a team might tank, it will tank when it decides it has no chance of making the playoffs
  # (after a minimum number of games have been played, currently set to half their games)
  # This is when it would not be enough for the team to win all its remaining games 
  # to have a win percentage at least as good as the last playoff team
  # (assuming that cutoff win percentage remains the same)
  #
  # mode = 1 or 2:
  # When two non-tanking teams or two tanking teams play each other, the better team wins with probability gamma
  # When a tanking team plays a non-tanking team, the tanking team always loses
  # 
  # mode = 3 or 4:
  # Variants of (Zermelo-)Bradley-Terry model used to determine who wins each game
  #
  # Assumptions:
  # 1. No simultaneous games
  # 2. Teams keep same relative true ranking throughout season
  # 3. No conferences / divisions
  # 4. No home/away games
  ###

  ## Set constants
  USE_MATH_ELIM = true # false = use effective elimination, true = use mathematical elimination
  CALC_MATH_ELIM = 2 # 0: do not calculate, 1: use heuristic only, 2: use MIP
  step_size = 1 / num_steps
  array_of_tanking_probabilities = 0:step_size:1
  num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
  num_games_total = num_rounds * num_games_per_round
  num_team_games = num_rounds * (num_teams - 1)
  max_games_remaining = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen

  breakpoint_game_for_draft = [round(breakpoint_list[i] * num_games_total) for i in 1:length(breakpoint_list)]

  ## Prepare output
  avg_kend = zeros(Float64, num_steps+1, length(breakpoint_list))
  #avg_kend_h2h = zeros(Float64, num_steps+1, length(breakpoint_list))
  avg_games_tanked = zeros(Float64, num_steps+1, length(breakpoint_list))
  avg_already_tank = zeros(Float64, num_steps+1, length(breakpoint_list))
  avg_eliminated = zeros(Float64, num_steps+1, num_games_total)
  avg_kend_gold = 0.0
  avg_kend_lenten = 0.0

  ## Set up for game order 
  #games = Array{Int64}(undef, num_rounds, num_games_per_round, 2) # games played in each round
  schedule = Matrix{Int64}(undef, num_rounds * num_games_per_round, 2) # schedule of games
  ord_games = Matrix{Int64}(undef, num_games_per_round, 2) # all the games played per round [internal]
  # Alternative to below is 
  # using Combinatorics
  # ord_games = collect(combinations(1:num_teams,2))
  curr_game_ind = 1 
  for i in 1:num_teams
    for j in (i+1):num_teams
      ord_games[curr_game_ind,1] = i
      ord_games[curr_game_ind,2] = j
      curr_game_ind = curr_game_ind + 1
    end
  end

  ## For ranking teams
  # stats[:,1] is team "name"
  # stats[:,2] is num wins
  # stats[:,3] is remaining
  # stats[:,4] is when team is eliminated (in terms of remaining games)
  # stats[:,5] is when team is mathematically eliminated (in terms of remaining games)
  # stats[:,6] is win percentage
  # stats[:,7] is indicator for whether team tanks
  # not implemented:  stats[:,7] is team rank
  size_of_stats = 7
  team_name_ind = 1
  num_wins_ind = 2
  games_left_ind = 3
  games_left_when_elim_ind = 4
  games_left_when_math_elim_ind = 5
  elim_ind = USE_MATH_ELIM ? games_left_when_math_elim_ind : games_left_when_elim_ind
  win_pct_ind = 6
  will_tank_ind = 7
  stats = Matrix{Any}(undef, num_teams, size_of_stats)
  draft_rank_of_team = Array{Any}(undef, num_teams, length(breakpoint_game_for_draft))
  #draft_rank_of_team_h2h = Array{Any}(undef, num_teams, length(breakpoint_game_for_draft))
  #draft_ranking = Array{Any}(undef, num_teams, size_of_stats, length(breakpoint_game_for_draft))
  #draft_ranking_row_index = Array{Any}(undef, num_teams, length(breakpoint_game_for_draft))

  ## Begin calculations
  for step_ind in 1:length(array_of_tanking_probabilities)
    tank_perc = array_of_tanking_probabilities[step_ind]
    print("Simulating season with $tank_perc ratio of teams tanking\n")
    for rep = 1:num_replications
      print("\tRepeat $rep/$num_replications, ratio $tank_perc\n")
      #print("\t\tavg_kend \t$(avg_kend[step_ind,:])\n")
      ##print("\t\tavg_kend_h2h \t$(avg_kend_h2h[step_ind,:])\n")
      ## Set up stats for current repeat
      num_eff_elim = 0
      num_math_elim = 0
      num_teams_tanking = 0
      num_games_tanked = 0
      outcome = zeros(Int, num_games_total)
      num_wins_since_elim = zeros(Int, num_teams)
      elimination_index = zeros(Int, num_teams) # game when was this team eliminated
      for i = 1:num_teams
        stats[i,team_name_ind] = i # team name
        stats[i,num_wins_ind] = 0 # num wins
        stats[i,games_left_ind] = num_team_games # num games left
        stats[i,games_left_when_elim_ind] = -1 # when team is eliminated (in terms of how many left)
        stats[i,games_left_when_math_elim_ind] = -1 # when team is mathematically eliminated (in terms of how many left)
        stats[i,win_pct_ind] = 0.0 # win percentage
        if rand() > tank_perc
          stats[i,will_tank_ind] = 0 # will not tank
        else
          stats[i,will_tank_ind] = 1 # will tank
        end
      end # set up stats array

      ## Set up initial ranking
      rank_of_team = sortperm(randn(num_teams)) # initial ranking (returns rank of team i)
      team_in_pos = Array{Int}(undef, num_teams) # inverse ranking (returns team that is in position i)
      #rank_of_team_h2h = sortperm(randn(num_teams)) #rank_of_team # initial ranking (returns rank of team i)
      #team_in_pos_h2h = Array{Int}(undef, num_teams) # inverse ranking (returns team that is in position i)
      for i = 1:num_teams
        team_in_pos[rank_of_team[i]] = i
        #team_in_pos_h2h[rank_of_team_h2h[i]] = i
      end
      h2h = []
      if return_h2h
        h2h = zeros(Int, num_teams, num_teams)
      end

      ## Set random game order for this repeat
      for round_ind = 1:num_rounds
        perm = sortperm(randn(num_games_per_round))
        start_round = num_games_per_round * (round_ind - 1) + 1
        end_round = num_games_per_round * round_ind
        schedule[start_round:end_round, :] = ord_games[perm,:]
      end

      ## Setup MIP
      best_outcomes = zeros(Int, num_teams, num_games_total)
      best_num_wins = zeros(Int, num_teams, num_teams)
      best_rank = zeros(Int, num_teams)
      num_mips = 0
      model = CALC_MATH_ELIM > 1 ? setupMIP(schedule, num_teams, num_playoff_teams, num_team_games, num_games_total) : 0

      ## Run one season
      for game_ind = 1:num_games_total
          # Find cutoff [do this every game - find last playoff team - set that as cutoff]
          # Tie-breaking is fewest games left
          last_team = team_in_pos[num_playoff_teams]
          cutoff_avg = stats[last_team,win_pct_ind]

          # Current teams playing
          i = schedule[game_ind, 1] #games[round_ind, round_game_ind, 1]
          j = schedule[game_ind, 2] #games[round_ind, round_game_ind, 2]

          # Decide who wins the game
          team_i_wins = teamWillWin(i, j, stats, gamma, true_strength, mode, games_left_ind, elim_ind, will_tank_ind)
          outcome[game_ind] = team_i_wins ? i : j
          if return_h2h
            h2h[i,j] = h2h[i,j] + team_i_wins
            h2h[j,i] = h2h[j,i] + !team_i_wins
          end

          # Check tanking
          team_i_is_tanking = teamIsTanking(i, stats, games_left_ind, elim_ind, will_tank_ind)
          team_j_is_tanking = teamIsTanking(j, stats, games_left_ind, elim_ind, will_tank_ind)
          if team_i_is_tanking || team_j_is_tanking
            num_games_tanked += 1
          end

          # Do updates
          for k in [i,j]
            team_k_wins = (k == i) ? team_i_wins : !team_i_wins
            stats[k,num_wins_ind] = stats[k,num_wins_ind] + team_k_wins
            stats[k,games_left_ind] = stats[k,games_left_ind] - 1 # one fewer game remaining
            stats[k,win_pct_ind] = stats[k,num_wins_ind] / (num_team_games - stats[k,games_left_ind]) # update current win pct
            rank_of_team, team_in_pos = updateRank(stats, rank_of_team, team_in_pos, k, team_k_wins, num_teams, win_pct_ind, games_left_ind, h2h) # update rank
            #rank_of_team_h2h, team_in_pos_h2h = updateRank(stats, rank_of_team_h2h, team_in_pos_h2h, k, team_k_wins, num_teams, win_pct_ind, games_left_ind, h2h) # update rank
            
            # If team k wins and has been eliminated
            if team_k_wins && teamIsEffectivelyEliminated(stats[k,num_wins_ind], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
              num_wins_since_elim[k] = num_wins_since_elim[k] + 1
            end
          end
          #print("($i,$j) Team $i wins? $team_i_wins\n")
          #display([1:num_teams team_in_pos rank_of_team stats[team_in_pos,win_pct_ind]])

          # When the breakpoint for choosing a playoff ranking has been reached, set the ranking
          for r = 1:length(breakpoint_game_for_draft) 
            if game_ind == breakpoint_game_for_draft[r]
              draft_rank_of_team[:,r] = rank_of_team
              #draft_rank_of_team_h2h[:,r] = rank_of_team_h2h
              avg_already_tank[step_ind, r] += num_teams_tanking / num_replications
              avg_games_tanked[step_ind, r] += num_games_tanked / num_replications
              #draft_ranking[:,:,r] = stats
              #draft_ranking_row_index[:,r] = row_index
            end
          end
          
          # Maybe team is eliminated after this round; again check critical game for i,j (game that team is eliminated)
          last_team = team_in_pos[num_playoff_teams]
          cutoff_avg = stats[last_team,win_pct_ind]
          for k in 1:num_teams
            if stats[k,games_left_when_elim_ind] >= 0 # check team is not already eliminated
              continue
            end
            is_eliminated = teamIsEffectivelyEliminated(stats[k,num_wins_ind], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
            if is_eliminated
              num_eff_elim += 1
              stats[k,games_left_when_elim_ind] = stats[k,games_left_ind]
              if (!USE_MATH_ELIM)
                num_teams_tanking += stats[k,will_tank_ind] == 1
                elimination_index[k] = game_ind
              end
            end
          end # set critical game for teams i and j

          # Check mathematical elimination
          if (CALC_MATH_ELIM > 0)
            updateHeuristicBestRank!(outcome[game_ind], game_ind, schedule, best_outcomes, best_num_wins, best_rank)
            if (CALC_MATH_ELIM > 1)
              fixOutcome!(model, outcome[game_ind], game_ind, schedule)
            end
            for k in 1:num_teams
              if stats[k,games_left_when_math_elim_ind] >= 0 # check team is not already eliminated
                continue
              end
              (is_eliminated, mips_used) = teamIsMathematicallyEliminated!(k, game_ind, 
                  schedule, stats, outcome, best_outcomes, best_num_wins, best_rank, model,
                  num_teams, num_playoff_teams, num_team_games, num_games_total,
                  CALC_MATH_ELIM, num_wins_ind, games_left_ind)
              num_mips += mips_used
              if is_eliminated
                num_math_elim += 1
                stats[k,games_left_when_math_elim_ind] = stats[k,games_left_ind]
                if (USE_MATH_ELIM)
                  num_teams_tanking += stats[k,will_tank_ind] == 1
                  elimination_index[k] = game_ind
                end
              end
            end
          end
          #print("Game $game_ind\tNum MIPs: $num_mips\tNum eff elim: $num_eff_elim\tNum math elim: $num_math_elim\n")

          ## Update elimination stats
          num_eliminated = USE_MATH_ELIM ? num_math_elim : num_eff_elim
          avg_eliminated[step_ind, game_ind] += num_eliminated / num_replications
      end # iterate over games
      ## end of a season

      ## Make sure teams that were eliminated actually did not make the playoffs
      for i = 1:num_teams
        if (best_rank[i] < 0)
          @assert(rank_of_team[i] > num_playoff_teams)
        end
      end

      ## Get non-playoff teams at end of season
      nonplayoff_teams = team_in_pos[num_playoff_teams+1:num_teams]
      #nonplayoff_teams_h2h = team_in_pos_h2h[num_playoff_teams+1:num_teams]
      for r = 1:length(breakpoint_game_for_draft)
        tmp_stats = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        tmp_stats[:,1] = nonplayoff_teams
        tmp_stats[:,2] = draft_rank_of_team[nonplayoff_teams, r]
        sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
        avg_kend[step_ind, r] += kendtau_sorted(sorted_ranking[:,1], true_strength, mode) / num_replications
        #avg_kend[step_ind, r] += kendtau(draft_ranking[nonplayoff_teams,:,r], win_pct_ind, true_strength, mode) / num_replications

        # Now repeat with h2h
        #tmp_stats = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        #tmp_stats[:,1] = nonplayoff_teams_h2h
        #tmp_stats[:,2] = draft_rank_of_team_h2h[nonplayoff_teams, r]
        #sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
        #avg_kend_h2h[step_ind, r] += kendtau_sorted(sorted_ranking[:,1], true_strength, mode) / num_replications
      end

      if (tank_perc == 0.0)
        ## Also compute Kendall tau distance for the Gold and Lenten methods
        ranking_gold = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        ranking_gold[:,1] = nonplayoff_teams
        ranking_gold[:,2] = -1 * num_wins_since_elim[nonplayoff_teams] # negative because teams with more wins need to be ranked worse (as they are given a _better_ draft pick)
        avg_kend_gold += kendtau(ranking_gold, 2, true_strength, mode) / num_replications

        ## For the Lenten ranking, we need to double check that the teams we said are eliminated did not make the playoffs
        ## Note that if a team is mathematically eliminated, then it is also effectively eliminated; the problem is the converse
        ## For this experiment, it is ``okay'' if we rank a team that is effectively eliminated before another
        ## if in reality it ends up being _mathematically_ eliminated after the other
        ## What do we do with the teams that do not make the playoffs, but were never effectively eliminated?
        ## We will rank them in reverse order as they stand at the end of the season
        tmp_elim_index = Matrix{Int}(undef, num_teams - num_playoff_teams, 2)
        tmp_elim_index[:,1] = nonplayoff_teams
        tmp_elim_index[:,2] = elimination_index[nonplayoff_teams] # not using negative, because we will sort high-to-low later
        for elim_ind = 1:num_teams - num_playoff_teams
          if tmp_elim_index[elim_ind,2] == 0 # was never eliminated
            # Team was not eliminated but did not make the playoffs
            # It is so far unranked from Lenten perspective
            # Should be ranked higher than teams eliminated earlier
            # Among the teams not eliminated, pretend that higher rank at end of season means it was eliminated later
            curr_team_ind = tmp_elim_index[elim_ind,1]
            tmp_elim_index[elim_ind,2] = num_games_total + 1 + (num_teams - rank_of_team[curr_team_ind])
          end
        end
        ranking_lenten = sortslices(tmp_elim_index, dims=1, by = x -> x[2], rev=true) # descending; having a higher elimination index means Lenten ranks the team higher (since it was eliminated later), i.e., it has a worse draft pick
        avg_kend_lenten += kendtau_sorted(ranking_lenten[:,1], true_strength, mode) / num_replications
      end
    end # do replications
  end # looping over tanking percentages

  #avg_kend_to_return = return_h2h ? avg_kend_h2h : avg_kend
  return avg_kend, avg_games_tanked, avg_already_tank, avg_eliminated, avg_kend_gold, avg_kend_lenten
end # simulate
