#################################
# (Zermelo-)Bradley-Terry model #
#################################
# Aleksandr M. Kazachkov
# Shai Vardi
###
import Distributions.Beta

"""
Look at number of wins team with rank i won against team of rank j in every year
to calculate MLE for Bradley-Terry model
"""
function BT_MLE(;data_dir="../data", num_teams=30)
  #years = ["games1314.csv", "games1415.csv", "games1516.csv", "games1617.csv", "games1718.csv", "games1819.csv"]
  years = ["games0405.csv", "games0506.csv", "games0607.csv", "games0708.csv", "games0809.csv", "games0910.csv", "games1011.csv", "games1213.csv", "games1314.csv", "games1415.csv", "games1516.csv", "games1617.csv", "games1718.csv", "games1819.csv"]
  #years = ["games0708.csv"]
  num_years = length(years)

  num_stats = 6 # name, wins, losses, games_left, elim, win_pct
  num_wins_ind = 2
  win_pct_ind = 6

  h2h_by_year = zeros(Float64, num_years, num_teams, num_teams)
  stats = zeros(Float64, num_years, num_teams, num_stats)

  ## Retrieve stats for each year
  # h2hYYYY will have number of times team i beat team j
  # so total games between teams i and j are given there too
  # We also collect stats1314 to identify the ranking of each team at the end of each season
  total_games_played = zeros(Int, num_teams, num_teams)
  h2h_final = zeros(Float64, num_teams, num_teams)
  win_pct_nba = zeros(Float64, num_teams, num_years)
  for yr = 1:num_years
    _, _, stats[yr, :, :], _, h2h_by_year[yr, :, :] = parseNBASeason(years[yr], breakpoint_list, data_dir)

    # Get ranking from this season
    team_in_pos = sortperm(stats[yr, :, win_pct_ind], rev=true) # inverse ranking (returns team that is in position i)
    rank_of_team = Array{Int}(undef, num_teams)
    for i = 1:num_teams
      rank_of_team[team_in_pos[i]] = i
    end

    # Update h2h totals
    games_played = zeros(Int, num_teams)
    for i = 1:num_teams
      team_i = team_in_pos[i]
      for j = i+1:num_teams
        team_j = team_in_pos[j]
        total_games_played[i,j] += h2h_by_year[yr, team_i, team_j] + h2h_by_year[yr, team_j, team_i]
        total_games_played[j,i] += h2h_by_year[yr, team_i, team_j] + h2h_by_year[yr, team_j, team_i]
        h2h_final[i,j] += h2h_by_year[yr, team_i, team_j]
        h2h_final[j,i] += h2h_by_year[yr, team_j, team_i]
        games_played[i] += total_games_played[i,j]
      end
    end

    # should be the case that sum(h2h_by_year[yr, rank_of_team, : ], dims=2) - stats[yr, rank_of_team, num_wins_ind] is all zeros
    #println("h2h_yr$yr = ", h2h_by_year[yr, :, :])
    #println("win_pct_yr$yr = ", stats[yr, rank_of_team, win_pct_ind])
    #println("rank_yr$yr = ", rank_of_team)
    #println("team_in_pos_yr$yr = ", team_in_pos)
    #println("num_wins_yr$yr = ", stats[yr, rank_of_team, num_wins_ind])

    # Get win pct for this year
    win_pct_nba[:, yr] = stats[yr, rank_of_team, win_pct_ind]
  end # loop over years
  win_pct_nba_avg = 100 * sum(win_pct_nba[1:num_teams,:], dims=2)[:,1] / num_years
  #println("win_pct_nba_avg = ", win_pct_nba_avg)
  #println("total_games_played = ", total_games_played)
  #println("h2h_final = ", h2h_final)

  # Normalize h2h_final
  # Also, calculate number of comparisons won by i
  W = zeros(Float64, num_teams) # this will be actually the sum of the win probabilities
  for i = 1:num_teams
    W[i] = sum([h2h_final[i,j] / (h2h_final[i,j] + h2h_final[j,i]) for j in 1:num_teams if j != i])
    #for j = i+1:num_teams
      #if total_games_played[i,j] > 0
      #  h2h_final[i,j] /= total_games_played[i,j]
      #  h2h_final[j,i] /= total_games_played[i,j]
      #else
      #  h2h_final[i,j] = 0
      #  h2h_final[j,i] = 0
      #end

      #if h2h_final[i,j] >= h2h_final[j,i]
      #  W[i] += 1
      #end
      #if h2h_final[i,j] <= h2h_final[j,i]
      #  W[j] += 1
      #end
    #end
  end
  #println("W = ", W) 
  #println("h2h_final = ", h2h_final)

  p = ones(Float64, num_teams) / num_teams
  eps = 1e-5
  step = 0
  num_steps = 1000
  while step < num_steps
    p_new = BT_MLE_step(h2h_final, W, p)

    # Calculate differences
    sum_diff = 0.0
    for i = 1:num_teams
      sum_diff += abs(p[i] - p_new[i])
    end
    if sum_diff < eps
      break
    end
    p = p_new / sum(p_new)
    step += 1
  end
  #println("Num steps: $step")

  p = sort(p, rev=true)
  #p = sort([sum(h2h_final[:,i])/(num_teams-1) for i in 1:num_teams], rev=true)
  return p / sum(p)
end # BT_MLE

function BT_MLE_step(h2h, W, p)
  num_teams = length(p)
  p_new = copy(p)
  for i = 1:num_teams
    val = 0
    for j = 1:num_teams
      if j != i
        val += (h2h[i,j] + h2h[j,i]) / (p[i] + p[j])
      end
    end
    p_new[i] = W[i] * 1/val
  end
  return p_new
end # BT_MLE_step

function BT_avg(;data_dir="../data")
  # Get win_pct_nba
  win_pct_nba = readdlm(string(data_dir, "/winpct.csv"), ',')
  num_header_rows = 1
  num_years_nba = size(win_pct_nba, 2)
  p = sum(win_pct_nba[num_header_rows+1:num_header_rows+num_teams,:], dims=2) / num_years_nba
  p = p[:,1]
  return p
end # BT_avg

function test_BT(;data_dir="../data", distr=Beta(2,5), strength=[])
  if (length(strength) != num_teams)
    num_repeats = 100000
    strength = sort(rand(distr, num_teams), rev=true)
    for i = 2:num_repeats
      strength += sort(rand(distr, num_teams), rev=true)
    end
    strength /= num_repeats
    strength /= strength[1]
  end
  println("strength = ", strength)

  pr = zeros(Float64, num_teams, num_teams)
  exp_pr = zeros(Float64, num_teams, num_teams)
  for i=1:num_teams
    for j=1:num_teams
      if i==j
        continue
      end
      pr[i,j] = strength[i] / (strength[i] + strength[j])
      exp_pr[i,j] = exp(strength[i]) / (exp(strength[i]) + exp(strength[j]))
    end
  end
  
  win_pct_est = 100*[sum(pr[i,:]) for i in 1:num_teams] / (num_teams-1)
  exp_win_pct_est = 100*[sum(exp_pr[i,:]) for i in 1:num_teams] / (num_teams-1)
  #println("exp_win_pct_est = ", 100*exp_win_pct_est)

  ## Data for comparison
  win_pct_nba = readdlm(string(data_dir, "/winpct.csv"), ',') # [year, team]
  num_header_rows = 1
  num_years = size(win_pct_nba, 2)
  win_pct_nba_avg = sum(win_pct_nba[num_header_rows+1:num_header_rows+num_teams,:], dims=2)[:,1] / num_years

  println("win_pct_est = ", win_pct_est)
  println("win_pct_nba = ", win_pct_nba_avg)

  ## Calculate MSE
  loss = 0
  num_pts = num_teams * num_years
  for pos = 1:num_teams
    curr_calc = win_pct_est[pos]
    for year = 1:num_years
      curr_real = win_pct_nba[num_header_rows + pos,year]
      loss += (curr_real-curr_calc)^2 / num_pts # mean squared error
    end # loop over years
  end # loop over teams
  println("loss: $loss")
end # test_BT
