####################
# Useful functions #
####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
#include("mathelim.jl")

## Ranking type
# 1: strict: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
# 2: ties: [1,5] \succ [6,10] \succ \cdots \succ [26,30]
# 3: BT_uniform: Bradley-Terry with P(i>j) = p_i / (p_i + p_j); must also set distribution, where default is each team gets a strength score from U[0,1]
# 4: BT_exponential: Bradley-Terry with P(i>j) = exp(p_i) / (exp(p_i) + exp(p_j)); must also set distribution, where default is each team gets a strength score from U[0,1]; can consider others such as, e.g., using Beta(alpha=2, beta=5)
@enum MODE_TYPES NONE=0 STRICT TIES BT_DISTR BT_EXPONENTIAL BT_ESTIMATED

function greaterThanVal(x, y, eps = 1e-7)
  return x > y + eps
end # greaterThanVal
function lessThanVal(x, y, eps = 1e-7)
  return x < y - eps
end # lessThanVal
function isVal(x, y, eps = 1e-7)
  return abs(x-y) < eps
end # isVal

function updateStats!(x, y, num_replications, avg_stat = 1, stddev_stat = 2, min_stat = 3, max_stat = 4)
  x[avg_stat] = updateMoment(1, x[avg_stat], y, num_replications)
  x[stddev_stat] = updateMoment(2, x[stddev_stat], y, num_replications)
  x[min_stat] = updateMin(x[min_stat], y)
  x[max_stat] = updateMax(x[max_stat], y)
end # updateStats
function updateMoment(m, x, y, n=1)
  return x += (y^m)/n
end # updateMoment
function updateMin(x, y, eps = 1e-7)
  if lessThanVal(y,x,eps)
    return y
  else
    return x
  end
end # updateMin
function updateMax(x, y, eps = 1e-7)
  if greaterThanVal(y,x,eps)
    return y
  else
    return x
  end
end # updateMax

"""
runDraftLottery: using the given odds, calculate the draft order

We assume that nonplayoff_teams are ranked in order of best to worst record

The function is recursive, and calls itself to calculate the next team by setting odds of selected team to 0
(flagged by -1 in the code),
and then reweighting odds 


---
* odds: length(nonplayoff_teams) array of doubles, typically in increasing order (last place has higher chance)

---
Return the draft order
"""
function runDraftLottery(nonplayoff_teams, odds, num_teams, num_teams_selected_by_lottery)
  if num_teams == 0
    return []
  end
  
  # Check if we should not use the draft lottery,
  # but rather return the teams in the reverse order
  if num_teams_selected_by_lottery <= 0
    draft_order = zeros(Int, num_teams)
    tmp_ind = 1
    for curr_team_ind = length(nonplayoff_teams):-1:1
      if lessThanVal(odds[curr_team_ind], 0)
        continue
      else
        draft_order[tmp_ind] = nonplayoff_teams[curr_team_ind]
        tmp_ind += 1
        if tmp_ind > num_teams
          break
        end
      end
    end
    return draft_order
  end
   
  normalization = sum(i->(i>=-1e-7 ? i : 0), odds)
  curr_rand = rand()
  curr_sum = 0.
  selected_team_ind = -1
  selected_team = -1
  for curr_team_ind = 1:length(nonplayoff_teams)
    if lessThanVal(odds[curr_team_ind], 0)
      continue
    end

    curr_sum += odds[curr_team_ind] / normalization
    if curr_sum >= curr_rand
      selected_team = nonplayoff_teams[curr_team_ind]
      selected_team_ind = curr_team_ind
      break
    end
  end

  odds[selected_team_ind] = -1
  rank_rest = runDraftLottery(nonplayoff_teams, odds, num_teams-1, num_teams_selected_by_lottery-1)
  draft_order = pushfirst!(rank_rest, selected_team)
  return draft_order
end # runDraftLottery

function teamAdvances(i, ranks, num_playoff_teams_per_conf, conf=[])
  ###
  # teamAdvances
  #   i: team index
  #   ranks: vector of team ranks
  #   num_playoff_teams_per_conf: number of playoff teams from each conference
  #   conf: which conference each team belongs to (if empty, then assumed one conference)
  #
  # Return true if team advances
  if length(conf) == 0
    return ranks[i] <= num_playoff_teams_per_conf
  else
    mask = map(k->k==conf[i], conf)
    conf_ranks = ranks[mask]
    conf_ranki = count(k->k<ranks[i], conf_ranks)
    return conf_ranki <= num_playoff_teams_per_conf
  end
end # teamAdvances

function teamIsTanking(i, stats, elim_ind=4, will_tank_ind=7)
	return stats[i,will_tank_ind] == 1 && stats[i,elim_ind] >= 0
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
	elseif mode == BT_DISTR
		return (team_i_str) / ((team_i_str) + (team_j_str))
	elseif mode == BT_EXPONENTIAL
		return exp(team_i_str) / (exp(team_i_str) + exp(team_j_str))
	elseif mode == BT_ESTIMATED
		return (team_i_str) / ((team_i_str) + (team_j_str))
	end
end # teamIsBetter

function teamWillWinNoTanking(i, j, gamma, true_strength, mode)
	# Figure out which team is better and update gamma correspondingly
  # When BT is used, teamIsBetter returns the probability that team i beats team j
	better_team = teamIsBetter(i, j, true_strength, mode)
	if mode == STRICT || mode == TIES
		# if better_team == 1, keep gamma as it is
		if better_team == 0
			# If they are the same ranking then they have an equal chance of winning
			gamma = 0.5
		elseif better_team == -1
			gamma = 1 - gamma
		end
	elseif mode == BT_DISTR || mode == BT_EXPONENTIAL || mode == BT_ESTIMATED
		gamma = better_team # this is the probabilty that team i beats team j in the BT model
	end
		
	# team i (< j) wins with probability gamma (if i is indeed better than j)
	if rand() < gamma
		return true
	else
		return false
	end
end # teamWillWinNoTanking

function teamWillWin(i, j, stats, gamma, true_strength=30:-1:1, mode=STRICT, elim_ind=4, will_tank_ind=7)
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
	# stats[:,4/5] is num games games left when team was effectively/mathematically eliminated
	# stats[:,6] is win percentage
	# stats[:,7] is indicator for whether team tanks
	###
	team_i_tanks = teamIsTanking(i, stats, elim_ind, will_tank_ind) 
	team_j_tanks = teamIsTanking(j, stats, elim_ind, will_tank_ind) 

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

"""
teamIsMathematicallyEliminated
Check whether team k is mathematically eliminated

Sets the best known rank each non-eliminated team can achieve
as well as the corresponding stats and outcomes

Parameters
---
  * k: team ind
  * t: game ind
  * conf: which conference each team belongs to (if empty, then assumed one conference)

Updates
  * model
  * best_outcomes
  * best_h2h
  * best_num_wins
  * best_rank

Returns whether team k is eliminated and if a MIP solve was used
"""
function teamIsMathematicallyEliminated!(k, t, schedule, stats, outcome, h2h,
    best_outcomes, best_h2h, best_num_wins, best_rank, model,
    num_teams, num_playoff_teams, num_team_games, num_games_total, 
    math_elim_mode = -2, wins_ind = 2, games_left_ind = 3, conf=[])
  mips_used = 0
  ## Return immediately if team k is not eliminated in the existing heuristic solution
  if best_rank[k] > 0 && best_rank[k] <= num_playoff_teams
    return false, mips_used
  end
  ## Also return when the team is already mathematically eliminated
  if best_rank[k] < 0
    return true, mips_used
  end
  ## Also check the num_playoff_teams == 0 case
  if num_playoff_teams == 0
    best_rank[k] = -1
    return true, mips_used
  end

  heur_outcome, heur_h2h, heur_num_wins, heur_rank = 
      heuristicBestRank(k, t, schedule, stats, outcome, h2h,
          best_outcomes, best_h2h, best_num_wins, best_rank,
          num_playoff_teams, wins_ind, games_left_ind)

  ## If the heuristic solution is enough, stop here; else, go on to the MIP
  if heur_rank[k] <= num_playoff_teams
    best_outcomes[k,:] = heur_outcome
    best_h2h[k,:,:] = heur_h2h
    best_num_wins[k,:] = heur_num_wins
    best_rank[k] = heur_rank[k]
  elseif math_elim_mode > 1
    #print("Game $t, Team $k: Running MIP. Best rank: ", best_rank[k], " Heur rank: ", heur_rank[k], "\n")
    W =  stats[k, wins_ind] + stats[k, games_left_ind]
    fixVariables!(model, k, t+1, W, schedule, math_elim_mode)
    #W, w, x, y, z = setIncumbent!(model, k, t, schedule, heur_outcome)
    if solveMIP!(model)
      checkMIP(model, num_playoff_teams, math_elim_mode)
      updateUsingMIPSolution!(model, k, t, schedule, h2h, W, num_playoff_teams,
          best_outcomes, best_h2h, best_num_wins, best_rank, math_elim_mode)
    else
      best_rank[k] = -1
    end
              ### START DEBUG
              if false
                ## Save the hard LP
                lp_file = MathOptFormat.LP.Model()
                MOI.copy_to(lp_file, backend(model))
                name = string("model", math_elim_mode, ".lp")
                MOI.write_to_file(lp_file, name)
              end
              ### END DEBUG
    resetMIP!(model, t+1, schedule, h2h, stats, math_elim_mode)
    mips_used = 1
  end

  ## Update other teams if possible (if initialized and not yet eliminated)
  if best_rank[k] > 0
    updateOthersUsingBestSolution!(k, t, schedule, num_playoff_teams,
        best_outcomes, best_h2h, best_num_wins, best_rank)
  end

  ## If team k has been mathematically eliminated, mark it so by putting best rank as -1
  if best_rank[k] > num_playoff_teams
    best_rank[k] = -1
  end
  return best_rank[k] == -1, mips_used
end # teamIsMathematicallyEliminated

function teamIsEffectivelyEliminated(num_wins, num_games_remaining, num_team_games, cutoff_avg, max_games_remaining)
	###
	# teamIsEffectivelyEliminated
	#
	# If current playoff cutoff avg remains, 
	# and team k wins all its remaining games, 
	# will the team make it into the playoffs?
	###
	if num_games_remaining < max_games_remaining # make sure enough games have been played
		if lessThanVal((num_wins + num_games_remaining) / num_team_games, cutoff_avg)
			return true
		end
	end
	return false
end # teamIsEffectivelyEliminated

"""
kendtau_sorted
"""
function kendtau_sorted(sorted_ranking, true_strength=30:-1:1, mode=STRICT, min_rank=1)
	num_teams = size(sorted_ranking,1)
	kt = 0
	for i = min_rank:num_teams
		for j = i+1:num_teams
			better_team = teamIsBetter(Int(sorted_ranking[i]), Int(sorted_ranking[j]), true_strength, mode)
			out_of_order = false
      if mode == STRICT || mode == TIES
				out_of_order = (better_team == -1)
      elseif mode == BT_DISTR || mode == BT_EXPONENTIAL || mode == BT_ESTIMATED
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

"""
kendtau
Computes Kendell tau (Kemeny) distance

Column 1 should be index of the team
First sort by win percentage (stats[:,win_pct_ind])
This yields the noisy ranking
We want to calculate the distance to the true ranking (1,...,n)
Kendell tau distance is number of unordered pairs {i,j} for which noisy ranking disagrees with true ranking
K(\tau_1, \tau_2) 
  = |{ (i,j) : i < j, 
      (\tau_1(i) < \tau_1(j) && \tau_2(i) > \tau_2(j))
        ||
      (\tau_1(i) > \tau_2(j) && \tau_2(i) < \tau_2(j)) }|
"""
function kendtau(stats, win_pct_ind = 5, true_strength = 30:-1:1, mode=STRICT, min_rank=1)
	num_teams = size(stats,1)
	num_stats = size(stats,2)
	tmp = randn(num_teams) # randomize order among teams that have same number of wins
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

"""
rankTeamsFromWinTotals

Parameters
---
  * `num_wins`
  * `be_optimistic`: (default true) assume teams can achieve best-possible ranking based on win total (if false, use worst-possible)

Returns list `rank_of_team`
"""
function rankTeamsFromWinTotals(num_wins, be_optimistic = true)
  num_teams = length(num_wins)
 
  ## Identify ranking of teams based on the outcomes
  sorted_teams = sortperm(num_wins, rev=true)
  rank_of_team = Array{Int}(undef, num_teams)
  if be_optimistic # use optimistic ordering
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
    end # loop over teams
  else # else, use pessimistic ordering
    rank_of_team[sorted_teams[num_teams]] = num_teams
    ct = 1
    for i = num_teams-1:-1:1
      last_team = sorted_teams[i+1]
      curr_team = sorted_teams[i]
      last_wins = num_wins[last_team]
      curr_wins = num_wins[curr_team]
      rank_of_team[curr_team] = rank_of_team[last_team]
      if last_wins < curr_wins
        rank_of_team[curr_team] -= ct
        ct = 1
      else
        ct += 1
      end
    end # loop over teams
  end
  return rank_of_team
end # rankTeamsFromWinTotals

"""
updateRank
Update strict ranking and inverse map

After updating win percentage for team i and j, find their new postions
Tie breaking is as follows:
  1. win pct
  2. if tied, then head-to-head record
  3. if tied, then fewest games left
  4. if tied, then flip an unbiased coin

Parameters
---
  * [in/out] rank_of_team: current strict ranking
  * [in/out] team_in_pos
  * [in] stats
  * [in] team_i
  * [in] team_i_wins
  * [in] num_teams
  * [in] win_pct_ind
  * [in] games_left_ind
  * [in] h2h

Return strict ranking of teams and inverse (team that is in each position)
"""
function updateRank(rank_of_team, team_in_pos, stats, team_i, team_i_wins, num_teams, win_pct_ind = 6, games_left_ind = 3, h2h = [])
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

"""
teamIsContender
Check whether team k is a contender (is neither eliminated from the playoffs, nor guaranteed to make the playoffs) after game t

*Currently only checks that team is not eliminated from the playoffs before time t*

Parameters
---
  * [in] k: team ind
  * [in] t: game ind
  * [in] schedule: 
  * [in] stats:
  * [in] outcome:
  * [in] h2h:
  * [in] math_elim_game_total: game in which team was mathematically eliminated (< 0: never eliminated)
  * [in] num_playoff_teams

Returns whether team k is a contender after game t
"""
function teamIsContender(k, t, schedule, stats, outcome, h2h, math_elim_game_total, num_playoff_teams, wins_ind = 2, games_left_ind = 3)
  is_contender = math_elim_game_total[k] < 0 || math_elim_game_total[k] > t # team k is not eliminated yet (by game t)
  if is_contender
    _, _, _, rank_of_team = heuristicWorstRank(k, t, schedule, stats, outcome, h2h, wins_ind, games_left_ind)
    is_contender = (rank_of_team[k] > num_playoff_teams)
  end
  return is_contender
end # teamIsContender
