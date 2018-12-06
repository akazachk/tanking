##################### 
# Tanking simulator #
#####################
# Aleksandr M. Kazachkov
# Shai Vardi
###
include("utility.jl")

function simulate(num_teams, num_teams_in_playoffs, num_rounds, num_repeats, num_steps, gamma, set_ranking, true_strength, mode)
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
	###

	## Set constants
	step_size = 1 / num_steps
	array_of_tanking_probabilities = 0:step_size:1
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
	num_games = num_rounds * num_games_per_round
	num_team_games = num_rounds * (num_teams - 1)
	max_games_remaining = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen

	cutoff_game_for_draft = [round(set_ranking[i] * num_games) for i in 1:length(set_ranking)]

	## Prepare output
	avg_kend = zeros(Float64, num_steps+1, length(set_ranking))
	avg_games_tanked = zeros(Float64, num_steps+1, length(set_ranking))
	avg_already_tank = zeros(Float64, num_steps+1, length(set_ranking))
	avg_eliminated = zeros(Float64, num_steps+1, num_games)
	avg_kend_gold = 0.0
	avg_kend_lenten = 0.0

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
	# stats[:,1] is team "name"
	# stats[:,2] is num wins
	# stats[:,3] is remaining
	# stats[:,4] is when team is eliminated(in terms of remaining games)
	# stats[:,5] is win percentage
	# stats[:,6] is indicator for whether team tanks
	#	stats[:,7] is team rank
	size_of_stats = 6
	games_left_ind = 3
	win_pct_ind = 5
	stats = Matrix{Any}(undef, num_teams, size_of_stats)
	draft_rank_of_team = Array{Any}(undef, num_teams, length(cutoff_game_for_draft))
	#draft_ranking = Array{Any}(undef, num_teams, size_of_stats, length(cutoff_game_for_draft))
	#draft_ranking_row_index = Array{Any}(undef, num_teams, length(cutoff_game_for_draft))

	## Begin calculations
	for step_ind in 1:length(array_of_tanking_probabilities)
		tank_perc = array_of_tanking_probabilities[step_ind]
		print("Simulating season with $tank_perc ratio of teams tanking\n")
		for rep = 1:num_repeats
			print("\tRepeat $rep/$num_repeats (ratio $tank_perc)\n")
			## Set up stats for current repeat
			num_eliminated = 0
			num_teams_tanking = 0
			num_games_tanked = 0
			for i = 1:num_teams
				stats[i,1] = i # team name
				stats[i,2] = 0 # num wins
				stats[i,games_left_ind] = num_team_games # num games left
				stats[i,4] = 0 # when team is eliminated (in terms of how many left)
				stats[i,win_pct_ind] = 0.0 # win percentage
				if rand() > tank_perc
					stats[i,6] = 0 # will not tank
				else
					stats[i,6] = 1 # will tank
				end
			end # set up stats array

			## Set up initial ranking
			rank_of_team = sortperm(randn(num_teams)) # initial ranking (returns rank of team i)
			team_in_pos = Array{Int}(undef, num_teams) # inverse ranking (returns team that is in position i)
			for i = 1:num_teams
				team_in_pos[rank_of_team[i]] = i
			end
			h2h = zeros(Int, num_teams, num_teams)
			num_wins_since_elim = zeros(Int, num_teams)
			elimination_index = zeros(Int, num_teams) # when was this team eliminated

			## Set random game order for this repeat
			for round_ind = 1:num_rounds
				perm = sortperm(randn(num_games_per_round))
				games[round_ind,:,:] = ord_games[perm,:]
			end

			## Run one season
			game_ind = 0
			for round_ind = 1:num_rounds
				for round_game_ind = 1:num_games_per_round
					game_ind += 1

					# Find cutoff [do this every game - find last playoff team - set that as cutoff]
					# Tie-breaking is fewest games left
					#stats, row_index = sortTeams(stats)
					#last_team = num_teams_in_playoffs
					#last_team = sorted_teams[num_teams_in_playoffs,1]
					last_team = team_in_pos[num_teams_in_playoffs]
					cutoff_avg = stats[last_team,win_pct_ind]

					# Current teams playing
					i = games[round_ind, round_game_ind, 1]
					j = games[round_ind, round_game_ind, 2]
					#row_i = row_index[games[round_ind, round_game_ind, 1]]
					#row_j = row_index[games[round_ind, round_game_ind, 2]]

					# Set critical game for i,j (game that team is eliminated)
					for k in [i,j]
						if stats[k,4] == 0 # check team has not already started tanking
							if teamIsEliminated(stats[k,2], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
								stats[k,4] = stats[k,games_left_ind]
								num_eliminated += 1
								num_teams_tanking += stats[k,6] == 1
								elimination_index[k] = num_eliminated
							end
						end
					end # set critical game for teams i and j

					# Decide who wins the game
					team_i_wins = teamWillWin(i, j, stats, gamma, true_strength, mode)
					h2h[i,j] = h2h[i,j] + team_i_wins
					h2h[j,i] = h2h[j,i] + !team_i_wins

					# Check tanking
					if teamIsTanking(i, stats) || teamIsTanking(j, stats) #(stats[i,6] * stats[i,4] + stats[j,6] * stats[j,4] > 0)
					  num_games_tanked += 1
					end

					# Do updates
					for k in [i,j]
						team_k_wins = (k == i) ? team_i_wins : !team_i_wins
          	stats[k,2] = stats[k,2] + team_k_wins
						stats[k,games_left_ind] = stats[k,games_left_ind] - 1 # one fewer game remaining
						stats[k,win_pct_ind] = stats[k,2] / (num_team_games - stats[k,games_left_ind]) # update current win pct
						rank_of_team, team_in_pos = updateRank(stats, rank_of_team, team_in_pos, k, team_k_wins, num_teams, win_pct_ind, games_left_ind, h2h) # update rank
						
						# If team k wins and has been eliminated
						if team_k_wins && teamIsEliminated(stats[k,2], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
							num_wins_since_elim[k] = num_wins_since_elim[k] + 1
						end
					end
					avg_eliminated[step_ind, game_ind] += num_eliminated / num_repeats
					#print("($i,$j) Team $i wins? $team_i_wins\n")
					#display([1:num_teams team_in_pos rank_of_team stats[team_in_pos,win_pct_ind]])

					# When the cutoff for choosing a playoff ranking has been reached, set the ranking
					for r = 1:length(cutoff_game_for_draft) 
						if game_ind == cutoff_game_for_draft[r]
							draft_rank_of_team[:,r] = rank_of_team
							avg_already_tank[step_ind, r] += num_teams_tanking / num_repeats
							avg_games_tanked[step_ind, r] += num_games_tanked / num_repeats
							#draft_ranking[:,:,r] = stats
							#draft_ranking_row_index[:,r] = row_index
						end
					end
					
					# Maybe team is eliminated after this round; again check critical game for i,j (game that team is eliminated)
					last_team = team_in_pos[num_teams_in_playoffs]
					cutoff_avg = stats[last_team,win_pct_ind]
					for k in [i,j]
						if stats[k,4] == 0 # check team has not already started tanking
							if teamIsEliminated(stats[k,2], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
								stats[k,4] = stats[k,games_left_ind]
								num_eliminated += 1
								num_teams_tanking += stats[k,6] == 1
								elimination_index[k] = num_eliminated
							end
						end
					end # set critical game for teams i and j
				end # iterate over num_games_per_round
			end # iterate over rounds
			## end of a season

			## Get non-playoff teams at end of season
			np_index = team_in_pos[num_teams_in_playoffs+1:num_teams]
			for r = 1:length(cutoff_game_for_draft)
				tmp_stats = Matrix{Int}(undef, num_teams - num_teams_in_playoffs, 2)
				tmp_stats[:,1] = np_index
				tmp_stats[:,2] = draft_rank_of_team[np_index, r]
				sorted_ranking = sortslices(tmp_stats, dims=1, by = x -> x[2], rev=false) # ascending, as already in order
				avg_kend[step_ind, r] += kendtau_sorted(sorted_ranking[:,1], true_strength, mode) / num_repeats
				#avg_kend[step_ind, r] += kendtau(draft_ranking[np_index,:,r], win_pct_ind, true_strength, mode) / num_repeats
			end

			## Also compute Kendall tau distance for the Gold and Lenten methods
			ranking_gold = Matrix{Int}(undef, num_teams - num_teams_in_playoffs, 2)
			ranking_gold[:,1] = np_index
			ranking_gold[:,2] = -1 * num_wins_since_elim[np_index] # negative because teams with more wins need to be ranked worse (as they are given a _higher_ draft pick)
			avg_kend_gold += kendtau(ranking_gold, 2, true_strength, mode) / num_repeats

			## For the Lenten ranking, we need to double check that the teams we said are eliminated did not make the playoffs
			## Note that if a team is mathematically eliminated, then it is also effectively eliminated; the problem is the converse
			## For this experiment, it is ``okay'' if we rank a team that is effectively eliminated before another
			## if in reality it ends up being _mathematically_ eliminated after the other
			## What do we do with the teams that do not make the playoffs, but were never effectively eliminated?
			## We will rank them in reverse order as they stand at the end of the season
			tmp_elim_index = Matrix{Int}(undef, num_teams - num_teams_in_playoffs, 2)
			tmp_elim_index[:,1] = np_index
			tmp_elim_index[:,2] = elimination_index[np_index] # not using negative, because we will sort high-to-low later
			for elim_ind = 1:num_teams - num_teams_in_playoffs
				if tmp_elim_index[elim_ind,2] == 0 # was never eliminated
					# Team was not eliminated but did not make the playoffs
					# It is so far unranked from Lenten perspective
					# Should be ranked higher than teams eliminated earlier
					# Among the teams not eliminated, pretend that higher rank at end of season means it was eliminated later
					curr_team_ind = tmp_elim_index[elim_ind,1]
					tmp_elim_index[elim_ind,2] = num_eliminated + 1 + (num_teams - rank_of_team[curr_team_ind])
				end
			end
			ranking_lenten = sortslices(tmp_elim_index, dims=1, by = x -> x[2], rev=true) # descending; having a higher elimination index means Lenten ranks the team higher (since it was eliminated later), i.e., it has a worse draft pick
			avg_kend_lenten += kendtau_sorted(ranking_lenten[:,1], true_strength, mode) / num_repeats
		end # do repeats
	end # looping over tanking percentages

	return avg_kend, avg_games_tanked, avg_already_tank, avg_eliminated, avg_kend_gold, avg_kend_lenten
end # simulate
