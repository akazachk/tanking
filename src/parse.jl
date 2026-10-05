###################################
# Parse NBA data from spreadsheet #
###################################
# Aleksandr M. Kazachkov
# Shai Vardi
###

#using CSV
#using DataFrames
using DelimitedFiles
#include("utility.jl")

## Eastern Conference franchises, under every name used since 2004-05; all other teams are in the West
## (e.g., "Charlotte Bobcats"/"Charlotte Hornets" are East, whereas "New Orleans Hornets" is West)
const EAST_TEAM_NAMES = Set([
  "Atlanta Hawks", "Boston Celtics", "Brooklyn Nets", "New Jersey Nets",
  "Charlotte Bobcats", "Charlotte Hornets", "Chicago Bulls", "Cleveland Cavaliers",
  "Detroit Pistons", "Indiana Pacers", "Miami Heat", "Milwaukee Bucks",
  "New York Knicks", "Orlando Magic", "Philadelphia 76ers", "Toronto Raptors",
  "Washington Wizards"
])

"""
    parseNBASeason

Replay an NBA regular season from `data_dir/filename`
(columns: Date, Time, Team1 = visitor, Score1, Team2 = home, Score2, ...),
tracking standings and effective elimination.

Canceled games (recorded with equal scores, e.g., BOS-IND on 2013-04-16) are skipped.
"""
function parseNBASeason(filename="games1314.csv", breakpoint_list=[3//4,1], data_dir=DATA_DIR)
	###
	# Parse real data
	###
	#df = CSV.read(filename, types = [String, String, String, Int, String, Int, String, Union{String,Missing}, Int, Union{String,Missing}])
	fullpath = joinpath(data_dir, filename)
	if !isfile(fullpath)
		error("NBA data file $fullpath not found; download it with `python3 scripts/fetch_bbref_games.py <season end year>`")
	end
	df = readdlm(fullpath, ',')
	num_header_rows = 1
	start_row = num_header_rows + 1
	# NB: the first team listed is the visitor and the second is the home team;
	# the names below are historical (home advantage is not used anywhere)
	home_team_ind = 3;
	away_team_ind = 5;
	home_score_ind = 4;
	away_score_ind = 6;

	num_rows = size(df)[1]
	num_cols = size(df)[2]
	num_games_total = num_rows	- num_header_rows
	cutoff_game_for_draft = [round(breakpoint_list[i] * num_games_total) for i in 1:length(breakpoint_list)]

	## Set constants
  USE_MATH_ELIM = false # false: use effective elimination, true: use mathematical elimination
  CALC_MATH_ELIM = 0 # 0: do not calculate, 1: use heuristic only, 2: use MIP
  h2h_left = []
  math_elim_mode = 0

	num_teams = 30
	num_teams_per_conf = 15
	num_playoff_teams_per_conf = 8

	num_team_games = 82 # scheduled games per team
	max_games_remaining = (4/5) * num_team_games # a minimum number of games needs to be played before tanking might happen
	# end constants

	## Teams and conferences (assigned by name, as alphabetical positions change when teams are renamed/relocated)
	teams = sort(unique(vcat(df[start_row:num_rows,home_team_ind], df[start_row:num_rows,away_team_ind])))
	@assert(length(teams) == num_teams, "Expected $num_teams teams in $filename, found $(length(teams))")
	is_east = [team in EAST_TEAM_NAMES for team in teams]
	east_teams = findall(is_east)
	west_teams = findall(.!is_east)
	@assert(length(east_teams) == num_teams_per_conf && length(west_teams) == num_teams_per_conf,
	        "Conference sizes are $(length(east_teams)) (East) and $(length(west_teams)) (West) in $filename; update EAST_TEAM_NAMES")

	## Canceled games are recorded with equal scores; count the games each team actually plays
	is_canceled = [df[row, home_score_ind] == df[row, away_score_ind] for row in start_row:num_rows]
	games_scheduled = zeros(Int, num_teams)
	for row = start_row:num_rows
		if !is_canceled[row - start_row + 1]
			games_scheduled[searchsortedfirst(teams, df[row, home_team_ind])] += 1
			games_scheduled[searchsortedfirst(teams, df[row, away_team_ind])] += 1
		end
	end

	## Prepare output
	num_eliminated = 0
	num_games_tanked = 0
	num_eliminated_by_game = zeros(Int, num_games_total)
	num_eliminated_at_cutoff = zeros(Int,length(breakpoint_list))
	num_games_tanked_at_cutoff = zeros(Int,length(breakpoint_list))

	## Set up data matrices
  schedule = df[start_row:num_rows,[home_team_ind,away_team_ind]]
  outcome = zeros(Int, num_games_total)
	#results = Matrix{Bool}(undef, num_teams, num_team_games)
	teamcounter = zeros(Int,num_teams,1) # number of games played by each team

  ## Set up MIP
  best_outcomes = zeros(Int, num_teams, num_games_total)
  best_num_wins = zeros(Int, num_teams, num_teams)
  best_rank = zeros(Int, num_teams)
  num_mips = 0
  model = CALC_MATH_ELIM > 1 ? setupMIP(schedule, h2h_left, num_teams, num_playoff_teams, num_team_games, num_games_total, math_elim_mode) : 0

	## Critical game is the point at which the team is eliminated (effectively or mathematically)
	## where effectively means that it cannot make the playoffs if the last team continues playing at the same average win rate
	critical_game = zeros(Int,num_teams,3) # [ elimination game, # wins @ elim, # losses @ elim ]

	## Set up stats
  size_of_stats = 6
	name_ind = 1
	wins_ind = 2
	losses_ind = 3
	games_left_ind = 4
  elim_ind = 5 # games left when eliminated
	win_pct_ind = 6
	stats = Matrix{Any}(undef, num_teams, size_of_stats) # [team name, wins, losses, games left, win pct]
	for i = 1:num_teams
		stats[i,name_ind] = i # name
		stats[i,wins_ind] = 0 # wins
		stats[i,losses_ind] = 0 # losses
		stats[i,games_left_ind] = games_scheduled[i] # games left
    stats[i,elim_ind] = -1 # when team is eliminated (in terms of how many left)
		stats[i,win_pct_ind] = 0.0 # win pct
	end

	## Set up initial ranking
	h2h = zeros(Int, num_teams, num_teams)
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

		## Check who wins the game
		home_score = df[row, home_score_ind]
		away_score = df[row, away_score_ind]
		canceled = is_canceled[game_ind]
		winner = (home_score > away_score) ? hometeam : awayteam
		loser = (home_score > away_score) ? awayteam : hometeam
		if !canceled
			outcome[game_ind] = winner
			h2h[winner,loser] += 1

			## Check if this game is tanked
			if (critical_game[hometeam,1] + critical_game[awayteam,1]) > 0
				num_games_tanked += 1
			end
		end

		## Set for each cutoff the number eliminated and number of games tanked
		for r = 1:length(breakpoint_list)
			if game_ind == cutoff_game_for_draft[r]
				num_eliminated_at_cutoff[r] = num_eliminated
				num_games_tanked_at_cutoff[r] = num_games_tanked
			end
		end
		num_eliminated_by_game[game_ind] = num_eliminated
		if canceled
			continue
		end

		## Do updates
		for k in [hometeam, awayteam]
			teamcounter[k] = teamcounter[k] + 1
			team_k_wins = (k == winner) ? true : false
			stats[k,wins_ind] += team_k_wins # one more win
			stats[k,losses_ind] += !team_k_wins # one more loss
			stats[k,games_left_ind] -= 1 # one fewer game remaining
			stats[k,win_pct_ind] = stats[k,wins_ind] / teamcounter[k] # update current win pct
			#results[k,teamcounter[k]] = team_k_wins
			if is_east[k]
				rank_of_team, team_in_pos_east = updateRank(rank_of_team, team_in_pos_east, stats, k, team_k_wins, num_teams_per_conf, win_pct_ind, games_left_ind, h2h) # update rank
			else
				rank_of_team, team_in_pos_west = updateRank(rank_of_team, team_in_pos_west, stats, k, team_k_wins, num_teams_per_conf, win_pct_ind, games_left_ind, h2h) # update rank
			end
		end

		## Check whether teams are eliminated
		#for k in [hometeam, awayteam]
    for k in 1:num_teams
			if stats[k,elim_ind] >= 0 ## check team is not already eliminated
        continue
      end
      if USE_MATH_ELIM
        # Mathematical elimination: solve MIP if heuristic does not find good schedule
        (is_eliminated, mips_used) = teamIsMathematicallyEliminated!(k, game_ind, 
            schedule, stats, outcome, best_outcomes, best_num_wins, best_rank, model,
            num_teams_per_conf, num_playoff_teams_per_conf, num_team_games, num_games_total,
            CALC_MATH_ELIM, wins_ind, games_left_ind)
         num_mips_used += mips_used
      else
        # Effective elimination: "If I win all my remaining games, and the cutoff for making the playoffs does not change, will I make the playoffs?"
        cutoff_avg = is_east[k] ? stats[team_in_pos_east[num_playoff_teams_per_conf],win_pct_ind] : stats[team_in_pos_west[num_playoff_teams_per_conf],win_pct_ind]
        is_eliminated = teamIsEffectivelyEliminated(stats[k,wins_ind], stats[k,games_left_ind], games_scheduled[k], cutoff_avg, max_games_remaining)
      end
      if is_eliminated
        stats[k,elim_ind] = stats[k,games_left_ind]
        critical_game[k,1] = teamcounter[k]
        critical_game[k,2] = stats[k,wins_ind]
        critical_game[k,3] = stats[k,losses_ind]
        num_eliminated += 1
      end
		end # check elimination

		#print("Game ",game_ind,": ",curr_home_name," (",hometeam,") vs ",curr_away_name," (",awayteam,"). Winner: ",winner,"\n")
	end # iterate over rows

	## Get non-playoff teams at the end of the season
	np_index_east = team_in_pos_east[num_playoff_teams_per_conf+1:num_teams_per_conf]
	np_index_west = team_in_pos_west[num_playoff_teams_per_conf+1:num_teams_per_conf]
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

	return num_eliminated_by_game, num_games_tanked_at_cutoff, stats, critical_game, h2h
end # parse
