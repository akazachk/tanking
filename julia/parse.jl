#using CSV
#using DataFrames
using DelimitedFiles
include("utility.jl")

function parseNBASeason(filename="games1314.xlsx", set_ranking=[3//4,1], data_dir="../data")
	###
	# Parse real data
	###
	#df = CSV.read(filename, types = [String, String, String, Int, String, Int, String, Union{String,Missing}, Int, Union{String,Missing}])
	df = readdlm(string(data_dir,"/",filename), ',')
	num_header_rows = 1
	start_row = num_header_rows + 1
	home_team_ind = 3;
	away_team_ind = 5;
	home_score_ind = 4;
	away_score_ind = 6;

	num_rows = size(df)[1]
	num_cols = size(df)[2]
	num_games = num_rows	- num_header_rows
	cutoff_game_for_draft = [round(set_ranking[i] * num_games) for i in 1:length(set_ranking)]

	## Set constants
	num_teams = 30
	num_teams_per_conf = 15
	num_teams_in_playoffs_per_conf = 8
	east_teams = [1 2 3 4 5 6 9 12 16 17 20 22 23 28 30]
	west_teams = [7 8 10 11 13 14 15 18 19 21 24 25 26 27 29]
	is_east = zeros(Bool, num_teams) # initialized to false
	for i in east_teams
		is_east[i] = true
	end

	num_team_games = 82 # ! set based on year
	max_games_remaining = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen
	# end constants

	## Prepare output
	num_eliminated = 0
	num_games_tanked = 0
	#num_eliminated_by_game = zeros(Int, num_games)
	num_eliminated_at_cutoff = zeros(Int,length(set_ranking))
	num_games_tanked_at_cutoff = zeros(Int,length(set_ranking))

	## Set up data matrices
	teams = sort(unique(df[start_row:num_rows,home_team_ind]))

	results = Matrix{Bool}(undef, num_teams, num_team_games)
	teamcounter = zeros(Int,num_teams,1) # number of games played by each team

	## Critical game is the point at which the team is eliminated,
	## in the sense that it cannot make the playoffs if the last team continues playing at the same average win rate
	critical_game = zeros(Int,num_teams,3) # [ crit game, # wins @ crit, # losses @ crit ]

	## Set up stats
	stats = Matrix{Any}(undef, num_teams, 5) # [team name, wins, losses, games left, win pct]
	name_ind = 1
	wins_ind = 2
	losses_ind = 3
	games_left_ind = 4
	win_pct_ind = 5
	for i = 1:num_teams
		stats[i,name_ind] = i # name
		stats[i,wins_ind] = 0 # wins
		stats[i,losses_ind] = 0 # losses
		stats[i,games_left_ind] = num_team_games # games left
		stats[i,win_pct_ind] = 0.0 # win pct
	end

	## Set up initial ranking
	rank_of_team = Array{Int}(undef, num_teams)
	rank_of_team[east_teams] = sortperm(randn(num_teams_per_conf))
	rank_of_team[west_teams] = sortperm(randn(num_teams_per_conf))
	team_in_pos_east = Array{Int}(undef, num_teams_per_conf) # inverse ranking (returns team that is in position i)
	team_in_pos_west = Array{Int}(undef, num_teams_per_conf) # inverse ranking (returns team that is in position i)
	for i = 1:num_teams_per_conf
		team_in_pos_east[rank_of_team[east_teams[i]]] = east_teams[i]
		team_in_pos_west[rank_of_team[west_teams[i]]] = west_teams[i]
	end

	## Iterate through the rows
	for row = start_row:num_rows
		game_ind = row - start_row + 1
		curr_home_name = df[row, home_team_ind]
		curr_away_name = df[row, away_team_ind]
		hometeam = searchsortedfirst(teams, curr_home_name) 
		awayteam = searchsortedfirst(teams, curr_away_name)

		## Decide who wins the game
		home_score = df[row, home_score_ind]
		away_score = df[row, away_score_ind]
		winning_team = (home_score > away_score) ? hometeam : awayteam
		losing_team = (home_score > away_score) ? awayteam : hometeam

		## Do updates
		for k in [hometeam, awayteam]
			teamcounter[k] = teamcounter[k] + 1
			team_k_wins = (k == winning_team) ? true : false
			stats[k,wins_ind] = stats[k,wins_ind] + team_k_wins # one more win
			stats[k,losses_ind] = stats[k,losses_ind] + !team_k_wins # one more loss
			stats[k,games_left_ind] = stats[k,games_left_ind] - 1 # one fewer game remaining
			stats[k,win_pct_ind] = stats[k,wins_ind] / teamcounter[k] # update current win pct
			results[k,teamcounter[k]] = team_k_wins
			if is_east[k]
				rank_of_team, team_in_pos_east = updateRank(stats, rank_of_team, team_in_pos_east, k, team_k_wins, num_teams_per_conf, win_pct_ind, games_left_ind) # update rank
			else
				rank_of_team, team_in_pos_west = updateRank(stats, rank_of_team, team_in_pos_west, k, team_k_wins, num_teams_per_conf, win_pct_ind, games_left_ind) # update rank
			end
		end

		## Critical game computation
		## "If I win all my remaining games, and the cutoff for making the playoffs does not change, will I make the playoffs?"
		cutoff_avg_home = 0
		cutoff_avg_away = 0
		if is_east[hometeam]
			cutoff_avg_home = stats[team_in_pos_east[num_teams_in_playoffs_per_conf],win_pct_ind]
		else
			cutoff_avg_home = stats[team_in_pos_west[num_teams_in_playoffs_per_conf],win_pct_ind]
		end
		if is_east[awayteam]
			cutoff_avg_away = stats[team_in_pos_east[num_teams_in_playoffs_per_conf],win_pct_ind]
		else
			cutoff_avg_away = stats[team_in_pos_west[num_teams_in_playoffs_per_conf],win_pct_ind]
		end

		## Set critical game for home and away teams
		for k in [hometeam, awayteam]
			if critical_game[k,1] == 0 ## check it has not been set yet
				cutoff_avg = (k == hometeam) ? cutoff_avg_home : cutoff_avg_away
				if teamIsEliminated(stats[k,wins_ind], stats[k,games_left_ind], num_team_games, cutoff_avg, max_games_remaining)
					critical_game[k,1] = teamcounter[k]
					critical_game[k,2] = stats[k,wins_ind]
					critical_game[k,3] = stats[k,losses_ind]
					num_eliminated += 1
				end
			end
		end

		## Check tanking
		if (critical_game[hometeam,1] + critical_game[awayteam,1]) > 0
			num_games_tanked += 1
		end

		## Set for each cutoff the number eliminated and number of games tanked
		for r = 1:length(set_ranking)
			if game_ind == cutoff_game_for_draft[r]
				num_eliminated_at_cutoff[r] = num_eliminated
				num_games_tanked_at_cutoff[r] = num_games_tanked
			end
		end
		#num_eliminated_by_game[game_ind] = num_eliminated

		#print("Game ",game_ind,": ",curr_home_name," (",hometeam,") vs ",curr_away_name," (",awayteam,"). Winner: ",winning_team,"\n")
	end # iterate over rows

	## Get non-playoff teams at the end of the season
	np_index_east = team_in_pos_east[num_teams_in_playoffs_per_conf+1:num_teams_per_conf]
	np_index_west = team_in_pos_west[num_teams_in_playoffs_per_conf+1:num_teams_per_conf]
	np_index = vcat(np_index_east, np_index_west)
	#print(critical_game[np_index_east,1])
	#print(critical_game[np_index_west,1],"\n")
	#for team_ind in np_index
	#	if critical_game[team_ind,1] == 0
	#		print(teams[team_ind]," was eliminated but was never tanking.\n")
	#	end
	#end

	#print("Teams not in the playoffs in the east and their win pct:\n\t", teams[np_index_east], "\n\t", stats[np_index_east,win_pct_ind],"\n"); 
	#print("Teams not in the playoffs in the west and their win pct:\n\t", teams[np_index_west], "\n\t", stats[np_index_west,win_pct_ind],"\n"); 

	return num_eliminated_at_cutoff, num_games_tanked_at_cutoff, stats, critical_game
end # parse
