##################### 
# Tanking simulator #
#####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
include("utility.jl")

function simulate(num_teams, num_rounds, num_repeats, num_steps, gamma, set_ranking)
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
	# When two non-tanking teams play each other, the better team wins with probability gamma
	# When a tanking team plays a non-tanking team, the tanking team always loses
	# When two tanking teams play each other, the one that is currently better wins
	#
	# We assume the "true" ranking is 1,...,n
	#
	# Assumptions:
	# 1. No simultaneous games
	# 2. Teams keep same relative true ranking throughout season
	# 3. No conferences
	# 4. Team decides whether to tank the next game after each game they play (not after other games)
	###

	## Set constants
	step_size = 1 / num_steps
	array_of_tanking_percentage = 0:step_size:1
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
	num_games = num_rounds * num_games_per_round
	num_team_games = num_rounds * (num_teams - 1)
	num_teams_in_playoffs = Int(2^ceil(log(2, num_teams / 2)))
	max_games_remaining = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen

	cutoff_game_for_draft = [round(set_ranking[i] * num_games) for i in 1:length(set_ranking)]
	avg_kend = Array{Float64}(undef, num_steps+1, length(set_ranking))
	avg_already_tank = Array{Float64}(undef, num_steps+1, length(set_ranking))

#	kend = Array{Int64}(undef, num_repeats, num_steps + 1, length(cutoff_game_for_draft))
# already_tank: number of tanking teams by the time draft ranking is set
#	already_tank = Array{Int64}(undef, num_repeats, num_steps + 1, length(cutoff_game_for_draft))

	## Set up for game order 
	games = Array{Int64}(undef, num_rounds, num_games_per_round, 2) # games played in each round
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
	# rank[:,1] is team "name"
	# rank[:,2] is num wins
	# rank[:,3] is remaining
	# rank[:,4] is critical game (in terms of remaining games)
	# rank[:,5] is win percentage
	# rank[:,6] is indicator for whether team tanks
	rank = Matrix{Any}(undef, num_teams, 6)
	draft_ranking = Array{Any}(undef, num_teams, 6, length(cutoff_game_for_draft))
	draft_ranking_row_index = Array{Any}(undef, num_teams, length(cutoff_game_for_draft))

	## Begin calculations
	for step_ind in 1:length(array_of_tanking_percentage)
		tank_perc = array_of_tanking_percentage[step_ind]
		print("Simulating season with $tank_perc ratio of teams tanking\n")
		for rep = 1:num_repeats
			print("\tRepeat $rep/$num_repeats\n")
			## Set up rank for current repeat
			for i = 1:num_teams
				rank[i,1] = i # team name
				rank[i,2] = 0 # num wins
				rank[i,3] = num_team_games # num games left
				rank[i,4] = 0 # critical game (in terms of how many left)
				rank[i,5] = 0.0 # win percentage
				if rand() > tank_perc
					rank[i,6] = 0 # will not tank
				else
					rank[i,6] = 1 # will tank
				end
			end # set up rank array

			## Set random game order for this repeat
			for round_ind = 1:num_rounds
				perm = sortperm(randn(num_games_per_round))
				#games[num_games_per_round * (round_ind - 1) + 1 : num_games_per_round * round_ind, :] = ord_games[perm,:]
				games[round_ind,:,:] = ord_games[perm,:]
			end

			## Run one season
			game_ind = 0
			for round_ind = 1:num_rounds
				for round_game_ind = 1:num_games_per_round
					game_ind += 1

					# Find cutoff [do this every game - find last playoff team - set that as cutoff]
					# Tie-breaking is fewest games left
					rank, row_index = sortTeams(rank)
					last_team = num_teams_in_playoffs
					#last_team = sorted_teams[num_teams_in_playoffs,1]
					cutoff_avg = rank[last_team, 5]

					# Current teams playing
					orig_i = games[round_ind, round_game_ind, 1]
					orig_j = games[round_ind, round_game_ind, 1]
					i = row_index[games[round_ind, round_game_ind, 1]]
					j = row_index[games[round_ind, round_game_ind, 2]]

					# Set critical game for i,j
					for k in [i,j]
						if (rank[k,6] == 1 && rank[k,4] == 0) # check team tanks and critical game is not set
							if setCriticalGame(rank[k,2], rank[k,3], 
									num_team_games, cutoff_avg, max_games_remaining)
								rank[k,4] = rank[k,3]
							end
						end
					end # set critical game for teams i and j

					# Decide who wins the game
					team_i_wins = teamWillWin(orig_i, orig_j, rank, gamma)
          rank[i,2] = rank[i,2] + team_i_wins
					rank[j,2] = rank[j,2] + !team_i_wins

					## Do updates
					# One fewer game remaining
					rank[i,3] = rank[i,3] - 1
					rank[j,3] = rank[j,3] - 1
					# Update current win pct
					rank[i,5] = rank[i,2] / (num_team_games - rank[i,3]) 
					rank[j,5] = rank[j,2] / (num_team_games - rank[j,3])

					# When the cutoff for choosing a playoff ranking has been reached, set the ranking
					for r = 1:length(cutoff_game_for_draft) 
						if game_ind == cutoff_game_for_draft[r]
							draft_ranking[:,:,r] = rank
							draft_ranking_row_index[:,r] = row_index
						end
					end
				end # iterate over num_games_per_round
			end # iterate over rounds
			## end of a season

			## Get non-playoff teams at end of season
			reverse_sorted, row_index = sortTeams(rank, false)
			np_index = reverse_sorted[1:num_teams-num_teams_in_playoffs]
			for r = 1:length(cutoff_game_for_draft)
				curr_row_index = draft_ranking_row_index[np_index,r]
				avg_kend[step_ind, r] += kendtau(draft_ranking[curr_row_index,:,r])
				avg_already_tank[step_ind, r] += sum([draft_ranking[i,4,r] > 0 for i in 1:num_teams])
			end
		end # do repeats
	
		## Compute average statistics over all repeats
		for r = 1:length(cutoff_game_for_draft)
			avg_kend[step_ind,r] /= num_repeats
			avg_already_tank[step_ind,r] /= num_repeats
		end
	end # looping over tanking percentages

	return avg_kend, avg_already_tank
end # simulate
