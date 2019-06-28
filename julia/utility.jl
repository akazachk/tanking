####################
# Useful functions #
####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
include("mathelim.jl")

## Ranking type
# 1: strict: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
# 2: ties: [1,5] \succ [6,10] \succ \cdots \succ [26,30]
# 3: BT_uniform: Bradley-Terry with P(i>j) = p_i / (p_i + p_j); must also set distribution, where default is each team gets a strength score from U[0,1]
# 4: BT_exponential: Bradley-Terry with P(i>j) = exp(p_i) / (exp(p_i) + exp(p_j)); must also set distribution, where default is each team gets a strength score from U[0,1]; can consider others such as, e.g., using Beta(alpha=2, beta=5)
@enum MODE_TYPES NONE=0 STRICT TIES BT_UNIFORM BT_EXPONENTIAL

function greaterThanVal(x, y, eps = 1e-7)
  return x > y + eps
end # greaterThanVal
function lessThanVal(x, y, eps = 1e-7)
  return x < y - eps
end # lessThanVal
function isVal(x, y, eps = 1e-7)
  return abs(x-y) < eps
end # isVal

function teamIsTanking(i, stats, games_left_ind=3, games_left_when_elim_ind=4, will_tank_ind=7)
	return stats[i,will_tank_ind] == 1 && stats[i,games_left_ind] <= stats[i,games_left_when_elim_ind]
end # teamIsTanking

function teamIsBetter(i, j, true_strength = 30:-1:1, mode=STRICT)
	###
	# teamIsBetter
	#
	# mode = 1: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
	# mode = 2: there may be ties in the teams, e.g., [1,5] \succ [6,10] \succ \cdots \succ [26,30]
	# mode = 3: use Bradley-Terry model, computing the probabilty team i beats team j via exponential score functions
	###
	team_i_str = true_strength[i]
	team_j_str = true_strength[j]
	if mode == STRICT || mode == TIES
		if team_i_str > team_j_str
			return 1
		elseif team_i_str < team_j_str
			return -1
		else
			return 0
		end
	elseif mode == BT_UNIFORM
		return (team_i_str) / ((team_i_str) + (team_j_str))
	elseif mode == BT_EXPONENTIAL
		return exp(team_i_str) / (exp(team_i_str) + exp(team_j_str))
	end
end # teamIsBetter

function teamWillWinNoTanking(i, j, gamma, true_strength, mode)
	# Figure out which team is better and update gamma correspondingly
	better_team = teamIsBetter(i, j, true_strength, mode)
	if mode == STRICT || mode == TIES
		# if better_team == 1, keep gamma as it is
		if better_team == 0
			# If they are the same ranking then they have an equal chance of winning
			gamma = 0.5
		elseif better_team == -1
			gamma = 1 - gamma
		end
	elseif mode == BT_UNIFORM || mode == BT_EXPONENTIAL
		gamma = better_team
	end
		
	# team i (< j) wins with probability gamma (if i is indeed better than j)
	if rand() < gamma
		return true
	else
		return false
	end
end # teamWillWinNoTanking

function teamWillWin(i, j, stats, gamma, true_strength=30:-1:1, mode=STRICT, games_left_ind=3, games_left_when_elim_ind=4, will_tank_ind=7)
	###
	# teamWillWin
	#
	# gamma is the probability the better team wins
	# mode = 1: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
	# mode = 2: there may be ties in the teams, e.g., [1,5] \succ [6,10] \succ \cdots \succ [26,30]
	# mode = 3: use Bradley-Terry model, computing the probabilty team i beats team j via exponential score functions
	#
	# Modes 1/2:
	# When two non-tanking teams play each other, the better team wins with probability gamma
	# When a tanking team plays a non-tanking team, the tanking team always loses; equivalent to setting gamma = 1
	# NB: a better model may be when the tanking team does so successfully with "some" probability
	# 2018/11/08: When two tanking teams play each other, it is treated as though neither is tanking
	# old way: When two tanking teams play each other, the one that is currently better wins
	#
	# Mode 3: 
	# Use Bradley-Terry model
	# 
	# stats[:,1] is team "name"
	# stats[:,2] is num wins
	# stats[:,3] is remaining
	# stats[:,4] is critical game (in terms of remaining games)
	# stats[:,5] is win percentage
	# stats[:,6] is indicator for whether team tanks
	###
	team_i_tanks = teamIsTanking(i, stats, games_left_ind, games_left_when_elim_ind, will_tank_ind) #stats[i,will_tank_ind] == 1 && stats[i,games_left_ind] <= stats[i,games_left_when_elim_ind] # team i is past the tanking cutoff point
	team_j_tanks = teamIsTanking(j, stats, games_left_ind, games_left_when_elim_ind, will_tank_ind) #stats[j,will_tank_ind] == 1 && stats[j,games_left_ind] <= stats[j,games_left_when_elim_ind] # team j is past the tanking cutoff point

	if (team_i_tanks && team_j_tanks) || (!team_i_tanks && !team_j_tanks)
		# Neither team is tanking, or both are; we treat this the same, as non-tanking
		return teamWillWinNoTanking(i, j, gamma, true_strength, mode)
	elseif team_i_tanks
		# Only team i tanks
		return false
	elseif team_j_tanks
		# Only team j tanks
		return true
	end # decide who wins the game
	return
end # teamWillWin

function teamIsMathematicallyEliminated!(k, t, schedule, stats, outcome,
    best_outcomes, best_num_wins, best_rank, model,
    num_teams, num_playoff_teams, num_team_games, num_games_total, 
    CALC_MATH_ELIM = 2, num_wins_ind = 2, games_left_ind = 3)
    ###
    # Check whether team k is mathematically eliminated
    #
    # Sets the best known rank each non-eliminated team can achieve
    # as well as the corresponding stats and outcomes
    #
    # Returns whether a MIP solve was used
    ###
    ## Return if team k is not eliminated in the existing heuristic solution
    mips_used = 0
    if best_rank[k] > 0 && best_rank[k] <= num_playoff_teams
      return false, mips_used
    end
    ## Also return when the team is already mathematically eliminated
    if best_rank[k] < 0
      return true, mips_used
    end

    heur_outcome, heur_num_wins, heur_rank = 
        heuristicBestRank(k, t, schedule, stats, outcome, num_wins_ind, games_left_ind)

    ## If the heuristic solution is enough, stop here; else, go on to the MIP
    if heur_rank[k] <= num_playoff_teams
      best_outcomes[k, :] = heur_outcome
      best_num_wins[k, :] = heur_num_wins
      best_rank[k] = heur_rank[k]
    elseif CALC_MATH_ELIM > 1
      #print("Game $t, Team $k: Running MIP. Best rank: ", best_rank[k], " Heur rank: ", heur_rank[k], "\n")
      fixVariables!(model, k, t+1, stats[k, num_wins_ind] + stats[k, games_left_ind], schedule)
      #W, w, x, y, z = setIncumbent!(model, k, t, schedule, heur_outcome)
      if solveMIP!(model, num_playoff_teams)
        updateUsingMIPSolution!(model, k, t, schedule, best_outcomes, best_num_wins, best_rank)
      else
        best_rank[k] = -1
      end
      resetMIP!(model, t+1, schedule, stats)
      mips_used = 1
    end

    ## Update other teams if possible
    if best_rank[k] > 0
      updateOthersUsingBestSolution!(k, t, schedule, best_outcomes, best_num_wins, best_rank)
    end

    ## If team k has been mathematically eliminated, mark it so by putting best rank as -1
    if best_rank[k] > num_playoff_teams
      best_rank[k] = -1
    end
    return best_rank[k] == -1, mips_used
end # teamIsMathematicallyEliminated

function teamIsEliminated(num_wins, num_games_remaining, num_team_games, cutoff_avg, max_games_remaining)
	###
	# teamIsEliminated
	#
	# If current playoff cutoff avg remains, 
	# and team k wins all its remaining games, 
	# will the team make it into the playoffs?
	###
	if num_games_remaining < max_games_remaining # make sure enough games have been played
		if (num_wins + num_games_remaining) / num_team_games < cutoff_avg - 1e-7
			return true
		end
	end
	return false
end # teamIsEliminated

function kendtau_sorted(sorted_ranking, true_strength=30:-1:1, mode=STRICT, min_rank=1)
	len = size(sorted_ranking,1)
	kt = 0
	for i = min_rank:len
		for j = i+1:len
			better_team = teamIsBetter(Int(sorted_ranking[i]), Int(sorted_ranking[j]), true_strength, mode)
			out_of_order = false
      if mode == STRICT || mode == TIES
				out_of_order = (better_team == -1)
      elseif mode == BT_UNIFORM || mode == BT_EXPONENTIAL
				out_of_order = (better_team < 0.5 - 1e-7)
			end
			#print("i: ", sorted_ranking[i], "\tj: ", sorted_ranking[j], "\tteamIsBetter: ",better_team,"\n")
			if out_of_order
				kt = kt + 1
			end
		end
	end
	return kt
end # kendtau_sorted

function kendtau(stats, win_pct_ind = 5, true_strength = 30:-1:1, mode=STRICT, min_rank=1)
	###
	# kendtau
	# Computes Kendell tau (Kemeny) distance
	#
	# Column 1 should be index of the team
	# First sort by win percentage (stats[:,win_pct_ind])
	# This yields the noisy ranking
	# We want to calculate the distance to the true ranking (1,...,n)
	# Kendell tau distance is number of unordered pairs {i,j} for which noisy ranking disagrees with true ranking
	# K(\tau_1, \tau_2) 
	#		= |{ (i,j) : i < j, 
	#				(\tau_1(i) < \tau_1(j) && \tau_2(i) > \tau_2(j))
	#					||
	#				(\tau_1(i) > \tau_2(j) && \tau_2(i) < \tau_2(j)) }|
	###
	len = size(stats,1)
	num_stats = size(stats,2)
	tmp = randn(len) # randomize order among teams that have same number of wins
	sorted_ranking = sortslices([stats tmp], dims=1, by = x -> (x[win_pct_ind],x[num_stats+1]), rev=true)
	return kendtau_sorted(sorted_ranking[:,1], true_strength, mode, min_rank)
end # kendtau

function sortTeams(stats, best_to_worst = true, win_pct_ind = 5, games_left_ind = 3)
	sorted = sortslices(stats, dims=1, by = x -> (x[win_pct_ind],-x[games_left_ind]), rev=best_to_worst) # sort descending (best-to-worst)
	#row_index = Array{Int}(undef, size(sorted)[1])
	#for i in 1:size(sorted)[1]
	#	row_index[sorted[i,1]] = i
	#end
	#return sorted, row_index
	return sorted
end # sortTeams

function updateRank(stats, rank_of_team, team_in_pos, team_i, team_i_wins, num_teams, win_pct_ind = 5, games_left_ind = 3, h2h = [])
	###
	# updateRank
	#
	# After updating win percentage for team i and j, find their new postions
	# Tie breaking is as follows:
	# 1. win pct
	# 2. if tied, then head-to-head record
	# 3. if tied, then fewest games left
	# 4. if tied, then flip an unbiased coin
	###
	old_rank_i = rank_of_team[team_i]
	win_pct = stats[team_i, win_pct_ind]
	games_left = stats[team_i, games_left_ind]
	inc = team_i_wins ? -1 : +1
	k = old_rank_i + inc
	should_continue = true
	while (k > 0 && k < num_teams+1 && should_continue)
		team_j = team_in_pos[k]
		curr_win_pct = stats[team_j, win_pct_ind]
		curr_games_left = stats[team_j, games_left_ind]
		ivj = (size(h2h,1) > 0) ? h2h[team_i,team_j] - h2h[team_j,team_i] : 0
		if team_i_wins
			should_continue = (win_pct > curr_win_pct) || 
												(win_pct == curr_win_pct && ivj > 0) ||
												(win_pct == curr_win_pct && ivj == 0 && games_left < curr_games_left) ||
												(win_pct == curr_win_pct && ivj == 0 && games_left == curr_games_left && rand() > 1/2)
		else
			should_continue = (win_pct < curr_win_pct) || 
												(win_pct == curr_win_pct && ivj < 0) ||
												(win_pct == curr_win_pct && ivj == 0 && games_left > curr_games_left) ||
												(win_pct == curr_win_pct && ivj == 0 && games_left == curr_games_left && rand() > 1/2)
		end
		if should_continue
			k += inc
		end
	end # find new ranking
	new_rank_i = k - inc # above loop stops at one above/below the new rank of i

	## Update all ranks that need to be updated
	tmp = team_in_pos[new_rank_i]
	team_in_pos[new_rank_i] = team_i
	rank_of_team[team_i] = new_rank_i
	for i = new_rank_i-inc:-inc:old_rank_i
		old = team_in_pos[i]
		team_in_pos[i] = tmp
		rank_of_team[tmp] = i
		tmp = old
	end
	return rank_of_team, team_in_pos
end # updateRank
