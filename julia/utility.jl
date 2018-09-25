####################
# Useful functions #
####################

function teamWillWin(i, j, rank, gamma=gamma)
	###
	# teamWillWin
	#
	# We assume that the true ranking is 1, ..., n
	#
	# Method 1:
	# When two non-tanking teams play each other, the better team wins with probability gamma
	# When a tanking team plays a non-tanking team, the tanking team always loses
	# When two tanking teams play each other, the one that is currently better wins
	# 
	# rank[:,1] is team "name"
	# rank[:,2] is num wins
	# rank[:,3] is remaining
	# rank[:,4] is critical game (in terms of remaining games)
	# rank[:,5] is win percentage
	# rank[:,6] is indicator for whether team tanks
	###
	team_i_tanks = rank[i,3] <= rank[i,4] # team i is past the tanking cutoff point
	team_j_tanks = rank[j,3] <= rank[j,4] # team j is past the tanking cutoff point
	if team_i_tanks && team_j_tanks
		# Both teams tank
		if rank[i,5] > rank[j,5] # team i is better
			return true
		else # team j is better
			return false
		end
	elseif team_i_tanks
		# Only team i tanks
		return false
	elseif team_j_tanks
		# Only team j tanks
		return true
	else
		# Neither team is tanking
		# Thus team i (< j) wins with probability gamma
		if rand() < gamma
			return true
		else
			return false
		end
	end # decide who wins the game
end # teamWillWin

function setCriticalGame(num_wins, num_games_remaining, num_team_games, cutoff_avg, max_games_remaining)
	###
	# setCriticalGame
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
end # setCriticalGame

function kendtau(rank)
	###
	# kendtau
	# Computes Kendell-Tau (Kemeny) distance
	#
	# First sort by win percentage (rank[:,5])
	# This yields the noisy ranking
	# We want to calculate the distance to the true ranking (1,...,n)
	# Kendell-Tau distance is number of unordered pairs {i,j} for which noisy ranking disagrees with true ranking
	# K(\tau_1, \tau_2) 
	#		= |{ (i,j) : i < j, 
	#				(\tau_1(i) < \tau_1(j) && \tau_2(i) > \tau_2(j))
	#					||
	#				(\tau_1(i) > \tau_2(j) && \tau_2(i) < \tau_2(j)) }|
	###
	len = size(rank,1)
	kt = 0
	tmp = randn(len) # randomize order among teams that have same number of wins
	noisy_rank = sortslices([rank tmp], dims=1, by = x -> (x[5],x[7]), rev=true)
	for m = 1:len
		for n = m+1:len
			if noisy_rank[m,1] > noisy_rank[n,1]
				kt = kt + 1
			end
		end
	end
	return kt
end # kendtau

function sortTeams(rank, best_to_worst = true, win_pct_ind = 5, games_left_ind = 3)
	sorted = sortslices(rank, dims=1, by = x -> (x[win_pct_ind],-x[games_left_ind]), rev=best_to_worst) # sort descending (best-to-worst)
	row_index = Array{Int}(undef, size(sorted)[1])
	for i in 1:size(sorted)[1]
		row_index[sorted[i,1]] = i
	end
	return sorted, row_index
end # sortTeams
