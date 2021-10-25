using Plots # has to come before Tanking due to a conflict with PyPlot and a fonts library
using Tanking # "using" adds symbols into namespace
import Gurobi  # for Gurobi.Env()
using JuMP
using Printf
using MathOptInterface
using Random
using Statistics
using DelimitedFiles
global GRB_ENV = Gurobi.Env()

#import Random # "import" keeps namespace clean
#Random.seed!(628) # for reproducibility; NOTE: should be set once somewhere before this function is run
num_teams = 30
num_playoff_teams = Int(2^ceil(log(2, num_teams / 2)))
num_players = 60
num_rounds = 3
num_steps = 1
gamma = 0.71425
breakpoint_list = [1.]
nba_odds_old = [.250, .199, .156, .119, .088, .063, .043, .028, .017, .011, .008, .007, .006, .005]
nba_odds_old = nba_odds_old[14:-1:1]
nba_odds_new = [.140, .140, .140, .125, .105, .090, .075, .060, .045, .030, .020, .015, .010, .005]
nba_odds_new = nba_odds_new[14:-1:1]
nba_odds_flat = [1. / 14. for i in 1:14]
nba_odds_list = [nba_odds_new, nba_odds_old, nba_odds_flat]
nba_num_lottery = [4, 3, 4]
math_elim_mode = 0
mode = Tanking.BT_DISTR
mode_list_name = ["BT"]
selected_steps = 1
initial_team_strength = [29.0,49.0,42.0,22.0,39.0,19.0,33.0,54.0,41.0,57.0,53.0,48.0,48.0,37.0,33.0,39.0,60.0,36.0,33.0,17.0,49.0,42.0,51.0,19.0,53.0,39.0,48.0,58.0,50.0,32.0]
initial_team_strength = sort!(initial_team_strength, rev = true)
results_dir = "."
csvext = ".csv"

function run_season(true_strength=nothing)
  if isnothing(true_strength)
    true_strength = num_teams:-1:1
  end
  num_replications = 1

  print("\n## Running one season. ##\nOutput is win_pct[step_ind, team_ind, stat], where step_ind is an integer corresponding to a repeat of the simulation with a different ratio of teams tanking, team_ind takes a value from 1 to num_teams, and stat is either avg_stat = 1, stddev_stat = 2, min_stat = 3, max_stat = 4.\n\n")
  win_pct = Tanking.simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, nba_num_lottery, true_strength, mode, math_elim_mode, selected_steps, GRB_ENV, true)

  ## Stats we keep
  avg_stat    = 1
  stddev_stat = 2
  min_stat    = 3
  max_stat    = 4
  num_stats   = 4
  prefix      = ["avg_", "stddev_", "min_", "max_"]
  return win_pct[1, :, avg_stat]
  display(win_pct[1, :, avg_stat])
  print("\n")
end # run_season

function simulate_current_system(num_replications = 1)
  # Loop over ten years, collect data like in run_season
  Random.seed!(628) # for reproducibility

  seasons = 10
  pick_Strength = zeros(Float64, num_players)
  std_team_wins = zeros((num_teams,seasons))
  season_avg_win_count = zeros((num_teams,seasons))
  avg_wins_per_season = zeros((num_teams,seasons))
  season_avg_win_count_sq = zeros((num_teams,seasons))
  for repl in 1:num_replications
    println("Replication: ", repl)
    updated_Team_Strength = copy(initial_team_strength)
    for k in 1:seasons
      win_pct_season = updated_Team_Strength/sum(updated_Team_Strength)
      strength_after_season = sort!(run_season(win_pct_season), rev=false)

      #nonplayoff_teams =fill(0.0, (1,14))
      #for i in 14:-1:1
      #  nonplayoff_teams[i]  = strength_after_season[i]
      #end
      nonplayoff_teams = sort(strength_after_season[1:num_teams-num_playoff_teams], rev=true)

      #println("\n## Strength after season ##")
      #println(strength_after_season)
      #println("\n## Nonnplayoff teams ##")
      #println(nonplayoff_teams)
      #println("\n## Draft odds ##")
      #println(nba_odds_new)

      draft_order = Tanking.runDraftLottery(nonplayoff_teams, copy(nba_odds_new), num_playoff_teams, 4)

      #println("\n## Draft order ##")
      #println(draft_order)
      #println("\n## Updated Team Strength 1 ##")
      #println(updated_Team_Strength)

      for j in 1:num_players
        pick_Strength[j] = 5*.974*exp(-0.05531*(j-1))
      end
      for i in 1:14
        updated_Team_Strength[i] =  draft_order[i]+pick_Strength[i]
      end
      for i in 1:14
        updated_Team_Strength[i] = updated_Team_Strength[i] + pick_Strength[(num_teams + i)]
      end
      for i in 15:num_teams
        updated_Team_Strength[i] =  updated_Team_Strength[i]+pick_Strength[i]+pick_Strength[(num_teams + i)]
      end

      updated_Team_Strength = 1230. * updated_Team_Strength / sum(updated_Team_Strength)
      updated_Team_Strength = sort!(updated_Team_Strength,rev = true)

      #println("\n## Updated Team Strength 2 ##")
      #println(updated_Team_Strength)

      println("\nSeason ", k)
      for i in 1:num_teams
        println(initial_team_strength[i]," -> ",updated_Team_Strength[i])
      end
      for i in 1:num_teams
        season_avg_win_count[i,k] = season_avg_win_count[i,k] +updated_Team_Strength[i]
      end
      for i in 1:num_teams
        season_avg_win_count_sq[i,k] = season_avg_win_count_sq[i,k] + (updated_Team_Strength[i])^2
      end
    end#end Season Loop
  end#End Repl Loop
  for k in 1:seasons
    println("Season: ", k)
    for i in 1:num_teams
      avg_wins_per_season[i,k] = season_avg_win_count[i,k]/(num_replications)
      std_team_wins[i,k] = (season_avg_win_count_sq[i,k]/num_replications) - avg_wins_per_season[i,k]^2
      std_team_wins[i,k] =sqrt(std_team_wins[i,k])
      println("Average wins team ", i, ": ", avg_wins_per_season[i,k])
      println("Standard Deviation of team ", i, ": ", std_team_wins[i,k])
    end
  end
  
  ## Comute Gini index
  println("\n## Gini coeff ##")
  gini_coeff = zeros(Float64, seasons)
  for k in 1:seasons
    gini_coeff[k] = gini(avg_wins_per_season[:,k])
    println(gini_coeff[k])
  end
  
  system = "_current_system"
  writedlm(string(results_dir, "/", "avg_", "wins_per_season", system, csvext), avg_wins_per_season, ',')
  writedlm(string(results_dir, "/", "avg_", "gini_coeff", system, csvext), gini_coeff, ',')
end #end simulate_current_system


function build_model(max_players::Int, player_value::Vector{Float64})
  model = Model(
                optimizer_with_attributes(() -> Gurobi.Optimizer(),
                                          "TimeLimit" => 10,
                                          "OutputFlag" => 1))

  @variable(model, x[i=1:num_teams,j=1:num_players], lower_bound = 0, upper_bound = 1,integer=true)

  b = fill(max_players, (1,num_teams)) # max_player limit
  
  for i in 1:num_teams
    con = @constraint(model, sum(x[i,j] for j=1:num_players)<= b[i])
  end
  for i in 1:num_teams
    con = @constraint(model, sum(x[i,j] for j=1:num_players) >= 1)
  end
  for i in 2:num_teams
    con = @constraint(model, sum(player_value[j]*x[i-1,j] for j=1:num_players)-sum(player_value[j]*x[i,j] for j=1:num_players) <= 0)
  end
  for j in 1:num_players
    con = @constraint(model, sum(x[i,j] for i=1:num_teams)<= 1)
  end

  return model, x
end # build_model


function simulate_opt_system(;num_replications::Int = 1, max_players = 3)
  ## Set up variables and constraints
  Random.seed!(628) # for reproducibility

  seasons = 10
  season_avg_win_count = zeros((num_teams,seasons))
  std_team_wins = zeros((num_teams,seasons))
  avg_wins_per_season = zeros((num_teams,seasons))
  season_avg_win_count_sq = zeros((num_teams,seasons))

  # Import draft class wins from draft distribution
  # Adjust for variability using montecarlo methods using rand function
  w = zeros(Float64, num_players) # array to hold pick strength
  for j in 1:num_players
    w[j] = 5*.974*exp(-0.05531*(j-1)) #+ rand((-1/j):(5/j))
  end
  #totalDraftWins = sum(w[j] for j=1:num_players)
  #println("Total Wins: ", totalDraftWins)

  model, x = build_model(max_players, w)
  for repl in 1:num_replications
    println("Replication: ", repl)
    previous_strength = copy(initial_team_strength)
    for k in 1:seasons
      win_pct_season = previous_strength/sum(previous_strength)

      s = run_season(win_pct_season)

      @objective(model,Max,sum(sum(w[j]*x[i,j]+s[i] for j=1:num_players) for i=1:num_teams))
      status=optimize!(model)
      #count = 0

      for i in 1:num_teams
        for j in 1:num_players
          if (JuMP.value.(x[i,j]) != 0)
            #println("x[",i,",",j,"] = ", JuMP.value.(x[i,j]))
            #println("w[",j,"] = ", JuMP.value.(w[j]))
            if JuMP.value.(x[i,j]) == 1
              s[i] = (s[i]) + (JuMP.value.(w[j]))
            end
          end
        end
      end
      sum_of_s = sum(s)
      println("\nSeason ", k)
      previous_strength = sort!(previous_strength,rev = true)
      s = sort!(s,rev = true)
      for i in 1:num_teams
        s[i] = (s[i]/sum_of_s)*1230
        println(previous_strength[i]," -> ",s[i])
      end

      for i in 1:num_teams
        previous_strength[i] = s[i]
      end
      #for i in 1:num_teams
      #rep_team_wins[i] = rep_team_wins[i]+s[i]
      #end
      for i in 1:num_teams
        season_avg_win_count[i,k] = season_avg_win_count[i,k] +s[i]
      end
      for i in 1:num_teams
        season_avg_win_count_sq[i,k] = season_avg_win_count_sq[i,k] + (s[i])^2
      end
    end# end season loop
  end# end repl loop
  for k in 1:seasons
    println("Season: ", k)
    for i in 1:num_teams
      avg_wins_per_season[i,k] = season_avg_win_count[i,k]/(num_replications)
      std_team_wins[i,k] = (season_avg_win_count_sq[i,k]/num_replications) - (avg_wins_per_season[i,k])^2
      println(std_team_wins[i,k])
      std_team_wins[i,k] =sqrt(std_team_wins[i,k])
      println("Average wins team ", i, ": ", avg_wins_per_season[i,k])
      println("Standard Deviation of team ", i, ": ", std_team_wins[i,k])
    end
  end
  
  ## Comute Gini index
  println("\n## Gini coeff ##")
  gini_coeff = zeros(Float64, seasons)
  for k in 1:seasons
    gini_coeff[k] = gini(avg_wins_per_season[:,k])
    println(gini_coeff[k])
  end
  
  system = string("_opt_rhs", max_players, "_system")
  writedlm(string(results_dir, "/", "avg_", "wins_per_season", system, csvext), avg_wins_per_season, ',')
  writedlm(string(results_dir, "/", "avg_", "gini_coeff", system, csvext), gini_coeff, ',')
end # end simulate_opt_system


function plot_wins_and_gini_index()
  # Read in data
  systems = ["current", "opt_rhs3", "opt_rhs2"]
  avg_wins_per_season = Vector{Matrix}()

  for sys_ind = 1:length(systems)
    system = string("_",systems[sys_ind],"_system")
    push!(avg_wins_per_season, readdlm(string(results_dir, "/", "avg_", "wins_per_season", system, csvext), ','))
    push!(gini_coeff, readdlm(string(results_dir, "/", "avg_", "gini_coeff", system, csvext), ','))
  end

  # This plots using GR
  gr() # Set the backend to GR
  ext = ".png"

  for sys_ind = 1:length(systems)
    p = Plots.scatter(transpose(avg_wins_per_season[sys_ind]),
                title = string("Average wins per season for each team (", systems[sys_ind], ")"),
                xlabel = "Season",
                xticks = 0:1:10,
                ylabel = "Average wins per season",
                legend = false)
    p = Plots.scatter!(0:0, transpose(initial_team_strength))
    Plots.savefig(string("AVGWinsPerSeasonGraph_", systems[sys_ind], ext)) # Saves the CURRENT_PLOT as a .png
  end

  # Next plot Gini index
  prilim_gini_coeff = gini(initial_team_strength)
  p = Plots.scatter([0], [prilim_gini_coeff], 
              title = string("Gini coefficient per season (", "current vs. optimization", ")"),
              xlabel = "Season",
              xticks = -1:1:10,
              ylabel = "Gini Coefficient",
              color = "black", 
              label = "Starting Gini Coefficient", 
              markershape = :rect,
              markersize = 5,
              legend=:topright)
  labels = ["Current system", "Opt system rhs 3", "Opt system rhs 2"]
  col   = ["red",   "orange", "green",   "blue",   "violet", "black", "gray"]
  style = ["solid", "dashed", "dashdot", "dotted", "dashed", "solid", "solid"]
  shape = [:+, :circle, :vline, :star4]
  shapesize = [5, 2, 5, 5, 5, 5, 5]
  for sys_ind = 1:2
    p = Plots.scatter!(gini_coeff[sys_ind], label = labels[sys_ind], color = col[sys_ind], markershape=shape[sys_ind], markersize=5)
  end
  Plots.savefig(string("GiniCoeff", ext)) # Saves the CURRENT_PLOT as a .png
end # plot_wins_and_gini_index


function gini(data::Vector{<:Real})
  n = length(data)
  total_sum = cumsum(sort(data))
  return (n + 1 - 2 * sum(total_sum) / total_sum[end]) / n
end


if abspath(PROGRAM_FILE) == @__FILE__
  if length(ARGS) < 1
    run_season()
    #elseif length(ARGS) == 1
    #  run_season(ARGS[1])
  else
    print("*** ERROR: Too many arguments provided.\n")
  end

  simulate_current_system()
  simulate_opt_system()
end
