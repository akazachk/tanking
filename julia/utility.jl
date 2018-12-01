####################
# Useful functions #
####################

function teamIsTanking(i, stats)
	return stats[i,6] == 1 && stats[i,3] <= stats[i,4]
end # teamIsTanking

function teamWillWin(i, j, stats, gamma=gamma)
	###
	# teamWillWin
	#
	# We assume that the true ranking is 1, ..., n
	#
	# Method 1:
	# When two non-tanking teams play each other, the better team wins with probability gamma
	# When a tanking team plays a non-tanking team, the tanking team always loses; equivalent to setting gamma = 1
	# NB: a better model may be when the tanking team does so successfully with "some" probability
	# 2018/11/08: When two tanking teams play each other, it is treated 
	# old way: When two tanking teams play each other, the one that is currently better wins
	#
	# Method 2 (not implemented):
	# Use Bradley-Terry model
	# 
	# stats[:,1] is team "name"
	# stats[:,2] is num wins
	# stats[:,3] is remaining
	# stats[:,4] is critical game (in terms of remaining games)
	# stats[:,5] is win percentage
	# stats[:,6] is indicator for whether team tanks
	###
	team_i_tanks = teamIsTanking(i, stats) #stats[i,6] == 1 && stats[i,3] <= stats[i,4] # team i is past the tanking cutoff point
	team_j_tanks = teamIsTanking(j, stats) #stats[j,6] == 1 && stats[j,3] <= stats[j,4] # team j is past the tanking cutoff point
	if (team_i_tanks && team_j_tanks) || (!team_i_tanks && !team_j_tanks)
		# Neither team is tanking, or both are; we treat this the same, as non-tanking
		# Thus team i (< j) wins with probability gamma
		if rand() < gamma
			return true
		else
			return false
		end
		# Both teams tank
		#if stats[i,5] > stats[j,5] # team i is better
		#	return true
		#else # team j is better
		#	return false
		#end
	elseif team_i_tanks
		# Only team i tanks
		return false
	elseif team_j_tanks
		# Only team j tanks
		return true
	#else
		# Neither team is tanking
		# Thus team i (< j) wins with probability gamma
		#if rand() < gamma
		#	return true
		#else
		#	return false
		#end
	end # decide who wins the game
	return
end # teamWillWin

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

function kendtau(stats, win_pct_ind = 5)
	###
	# kendtau
	# Computes Kendell tau (Kemeny) distance
	#
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
	kt = 0
	tmp = randn(len) # randomize order among teams that have same number of wins
	noisy_stats = sortslices([stats tmp], dims=1, by = x -> (x[win_pct_ind],x[num_stats+1]), rev=true)
	for m = 1:len
		for n = m+1:len
			if noisy_stats[m,1] > noisy_stats[n,1]
				kt = kt + 1
			end
		end
	end
	return kt
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

function updateRank(stats, rank_of_team, team_in_pos, team_i, team_i_wins, num_teams, win_pct_ind = 5, games_left_ind = 3)
	###
	# updateRank
	#
	# After updating win percentage for team i and j, find their new postions
	###
	old_rank_i = rank_of_team[team_i]
	win_pct = stats[team_i, win_pct_ind]
	games_left = stats[team_i, games_left_ind]
	inc = team_i_wins ? -1 : +1
	k = old_rank_i + inc
	should_continue = true
	while (k > 0 && k < num_teams+1 && should_continue)
		curr_win_pct = stats[team_in_pos[k], win_pct_ind]
		curr_games_left = stats[team_in_pos[k], games_left_ind]
		if team_i_wins
			should_continue = (win_pct > curr_win_pct) || (win_pct == curr_win_pct && games_left < curr_games_left)
		else
			should_continue = (win_pct < curr_win_pct) || (win_pct == curr_win_pct && games_left > curr_games_left)
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
