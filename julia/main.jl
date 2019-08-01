######################################
# On Tanking and Competitive Balance #
######################################
# Aleksandr M. Kazachkov
# Shai Vardi
###
# June 2019
###
# Implemented are the bilevel, Gold, Lenten, and NBA lottery-based rankings (both pre- and post-2019 changes)
# 1. Bilevel: set a breakpoint; rank nonplayoff teams based on their relative order at the breakpoint
# 2. Gold: rank based on number of wins after elimination
# 3. Lenten: reverse of order in which teams were eliminated
# 4. NBA lottery-based system: see odds below, which are applied to the end-of-season standings
###

## Required dependencies
using DelimitedFiles
using LaTeXStrings
using Printf
include("simulate.jl")
include("parse.jl")

## NBA draft odds (for non-playoff teams, in reverse order)
# If teams are tied, then those teams receive odds that are the average of the odds for the positions they occupy
# Keep these in order that teams are ranked (not for the draft, but by wins)
nba_odds_old = [.250, .199, .156, .119, .088, .063, .043, .028, .017, .011, .008, .007, .006, .005]
nba_odds_new = [.140, .140, .140, .125, .105, .090, .075, .060, .045, .030, .020, .015, .010, .005]
nba_odds_flat = [1. / 14. for i in 1:14]
nba_odds_old = nba_odds_old[14:-1:1]
nba_odds_new = nba_odds_new[14:-1:1]
nba_odds_list = [nba_odds_new, nba_odds_old, nba_odds_flat]

## When the (draft) ranking will be set as fraction of number games
#breakpoint_list = [4//8; 5//8; 6//8; 7//8; 1]
breakpoint_list = [1//2; 2//3; 3//4; 5//6; 7//8; 1]
num_rankings = length(breakpoint_list)
shape = [:vline, :utriangle, :rect, :x, :triangle, :circle]
col = ["red", "orange", "green", "blue", "violet", "black"]
#color_for_cutoff_point = ["c" "b" "m" "r" "k"]
num_teams = 30 # number of teams
num_playoff_teams = Int(2^ceil(log(2, num_teams / 2)))

## Set ranking type
# 1: strict: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
# 2: ties: [1,5] \succ [6,10] \succ \cdots \succ [26,30]
# 3: BT_uniform: Bradley-Terry with P(i>j) = p_i / (p_i + p_j); must also set distribution, where default is each team gets a strength score from U[0,1]
# 4: BT_exponential: Bradley-Terry with P(i>j) = exp(p_i) / (exp(p_i) + exp(p_j)); must also set distribution, where default is each team gets a strength score from U[0,1]; can consider others such as, e.g., using Beta(alpha=2, beta=5)
MODE = STRICT

ranking_type = ""
true_strength = []
csvext = ".csv"
ext_folder = "pdf"
lowext_folder = "png"
ext = string(ranking_type,".",ext_folder)
lowext = string(ranking_type,"_low",".",lowext_folder)

function set_mode(mode=MODE)
  ### Set global variables based on which mode is being used
	if mode == NONE
		cssvext = ".csv"
	elseif mode == STRICT
		# 1 \succ 2 \succ \cdots \succ 30
		global ranking_type="_strict"
		global true_strength = num_teams:-1:1
	elseif mode == TIES
		# [1,5] \succ [6,10] \succ \cdots \succ [26,30]
		global ranking_type="_ties"
		global true_strength = [Int(ceil(i/5)) for i in num_teams:-1:1] # allows for ties # old: [i:i+4 for i in 1:5:num_teams-4]
	elseif mode == BT_UNIFORM
		# Options to consider:
		# uniform distribution (same as beta(1,1))
		# --> essentially perfectly imbalanced
		# nonuniform distribution, beta(2,2), or maybe we should do some kind of bimodal distribution
		# --> more weight on middle teams, smaller probability of very weak or very strong teams
		global ranking_type="_BT_uniform"
		global true_strength = rand(30,1)
	elseif mode == BT_EXPONENTIAL
		global ranking_type="_BT_exponential"
		global true_strength = rand(30,1)
	end
	global csvext = string(ranking_type,".csv")
	global ext = string(ranking_type,".",ext_folder)
	global lowext = string(ranking_type,"_low",".",lowext_folder)
end # set_mode

## For plotting
DO_PLOTTING=true
environment = read(`uname`, String)
if chomp(environment) != "Darwin"
	DO_PLOTTING=false
end
TITLE_FONTSIZE=10
AXIS_TITLE_FONTSIZE=10
TICK_LABEL_FONTSIZE=8
LEGEND_FONTSIZE=8
LEGEND_TITLE_FONTSIZE=8
DPI=200
use_pyplot = true
upscale = 1 # upscaling in resolution
if !use_pyplot
	#ext = ".svg"
	using Plots
	using StatsPlots
	#gr(dpi=DPI); # Pkg.add("GR")
	pyplot(dpi=DPI)
	#pgfplots(dpi=DPI)

	# Set defaults
	default(tick_direction=:out)
	default(titlefont=TITLE_FONTSIZE)
	default(tickfont=TICK_LABEL_FONTSIZE)
	default(legendfont=LEGEND_FONTSIZE)
	default(grid=false)
	#default(size=(600*upscale,400*upscale)) # plot canvas size
	#fntsm = Plots.font("sans-serif", 8.0 * upscale)
	#fntlg = Plots.font("sans-serif", 12.0 * upscale)
	#default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
else
	if DO_PLOTTING
		using PyCall
		pygui(:qt5) # others do not work on mac
		#PyCall.PyDict(matplotlib["rcParams"])["font.serif"] = ["Cambria"]
		using PyPlot
		rc("text", usetex=true)
		rc("font", family="serif")
		rc("axes.spines", right=false, top=false)
		rc("axes", titlesize=TITLE_FONTSIZE * upscale)
		rc("axes", labelsize=AXIS_TITLE_FONTSIZE * upscale)
		rc("xtick", labelsize=TICK_LABEL_FONTSIZE * upscale)
		rc("ytick", labelsize=TICK_LABEL_FONTSIZE * upscale)
		rc("legend", fontsize=LEGEND_FONTSIZE * upscale)
		rc("legend", title_fontsize=LEGEND_TITLE_FONTSIZE * upscale)
		rc("legend", labelspacing=0.25)
		rc("lines", linewidth=1 * upscale)
		rc("lines", solid_capstyle="round")
		#rc("figure", figsize=[6.4,4.8] * upscale)
		rc("figure", figsize=[6*upscale,4*upscale]) # note that axes may change depending on label size
		#rc("figure", figsize=[6*1.5,4*1.5])
		rc("savefig", transparent=false)
		rc("savefig", bbox="tight")
		rc("savefig", pad_inches=0.0015 * upscale) # to allow for g,y,f to be not cut off
		rc("savefig", dpi=DPI)
	end
end

function main_simulate(;do_simulation = true, num_replications = 100000, 
    do_plotting=true, mode=MODE, results_dir = "../results", 
    num_rounds = 3,  num_steps = 20, gamma = 0.75, 
    math_elim_mode = 2)
  ###
  # main_simulate: Simulate a season and plot output
  #   * do_simulation: when false, read data from files in results_dir
  #   * num_replications: how many times to simulate each data point
  #   * do_plotting: if false, only gather data, without plotting it
  #   * mode: which kind of true ranking is used
  #     1 or 2: 
  #       When two non-tanking teams or two tanking teams play each other, 
  #       the better team wins with probability gamma
  #       When a tanking team plays a non-tanking team, the tanking team always loses
  #     3 or 4:
  #       Variants of (Zermelo-)Bradley-Terry model used to determine who wins each game
  #   * results_dir: where results should be saved and can be found
	#   * num_rounds: a round consists of each team playing each other team
  #   * num_steps: discretization of [0,1] for tanking probability
	#   * gamma: probability a better-ranked team wins over a worse-ranked team
  #   * math_elim_mode: 
  #     0: use effective elimination, 
  #     1: use mathematical elimination, but calculated by heuristics only, 
  #     2: use math elim, binary MIP, team-wise formulation
  #     3: use math elim, general integer MIP, team-wise formulation
  #     4: use math elim, binary MIP, cutoff formulation
  #     5: use math elim, general integer MIP, cutoff formulation
  #     <0: use effective elimination, but calculate mathematical elimination
  ###
	set_mode(mode)

	## Variables that need to be set
	## end variables that need to be set

	## Set constants
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
	num_games = num_rounds * num_games_per_round

	## For output
	kend = 0 # [step,cutoff,stat], holds KT distance for each tanking probability and cutoff for draft ranking
	games_tanked = 0 # [step,cutoff,stat], number games tanked by cutoff
	already_tank = 0 # [step,cutoff,stat], number teams already tanking by the cutoff
	math_eliminated = 0 # [step,game,stat], number teams eliminated by each game
	eff_eliminated = 0 # [step,game,stat], number teams eliminated by each game
  kend_nba = 0 # [step,odds,stat] KT distance based on NBA ranking
	kend_gold = 0 # [:,stat] ranking based on number of wins since elimination point
	kend_lenten = 0 # [:,stat], ranking in order of first-to-eliminated

  ## Stats we keep
  avg_stat = 1
  stddev_stat = 2
  min_stat = 3
  max_stat = 4
  num_stats = 4
  prefix = ["avg_", "stddev_", "min_", "max_"]

	## Do simulation or retrieve data
	if do_simulation
		## Do simulation
    kend, games_tanked, already_tank, math_eliminated, eff_eliminated, kend_nba, kend_gold, kend_lenten = simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, true_strength, mode, math_elim_mode)

    for stat = 1:num_stats
      writedlm(string(results_dir, "/", prefix[stat], "kend", csvext), kend[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "games_tanked", csvext), games_tanked[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "already_tank", csvext), already_tank[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "math_eliminated", csvext), math_eliminated[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "eff_eliminated", csvext), eff_eliminated[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "kend_nba", csvext), kend_nba[:,:,stat], ',')
      #writedlm(string(results_dir, "/", prefix[stat], "kend_gold", csvext), kend_gold[:,stat], ',')
      #writedlm(string(results_dir, "/", prefix[stat], "kend_lenten", csvext), kend_lenten[:,stat], ',')
    end
    writedlm(string(results_dir, "/", "kend_gold", csvext), kend_gold, ',')
    writedlm(string(results_dir, "/", "kend_lenten", csvext), kend_lenten, ',')
	else
    ## Resize things
    num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
    num_games_total = num_rounds * num_games_per_round
    kend = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    games_tanked = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    already_tank = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    math_eliminated = zeros(Float64, num_steps+1, num_games_total, num_stats)
    eff_eliminated = zeros(Float64, num_steps+1, num_games_total, num_stats)
    kend_nba = zeros(Float64, num_steps+1, length(nba_odds_list), num_stats)
    kend_gold = zeros(Float64, 1, num_stats)
    kend_lenten = zeros(Float64, 1, num_stats)
    for stat = 1:num_stats
      kend[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "kend", csvext), ',')
      games_tanked[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "games_tanked", csvext), ',')
      already_tank[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "already_tank", csvext), ',')
      math_eliminated[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "math_eliminated", csvext), ',')
      eff_eliminated[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "eff_eliminated", csvext), ',')
      kend_nba[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "kend_nba", csvext), ',')
      #kend_gold[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "kend_gold", csvext), ',')
      #kend_lenten[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "kend_lenten", csvext), ',')
    end
    kend_gold = readdlm(string(results_dir, "/", "kend_gold", csvext), ',')
    kend_lenten = readdlm(string(results_dir, "/", "kend_lenten", csvext), ',')
		num_steps = size(kend,1) - 1
	end

	if (do_plotting)
		## Plot avg_kend (Kendell tau distance) for the bilevel ranking
		print("Plotting avg_kend: average swap distance\n")
		minx = 0
		incx = 0.1
		maxx = 1
		miny = Int(floor(findmin(kend[:,:,1])[1]))
		incy = 1
		maxy = Int(ceil(findmax(kend[:,:,1])[1]))
		titlestring = L"\mbox{Effect of $\delta$ on bilevel ranking of non-playoff teams}"
		xlabelstring = L"\mbox{Probability of tanking once eliminated}"
		ylabelstring = L"\mbox{Distance from true ranking of non-playoff teams}"
		legendtitlestring = L"\mbox{Breakpoint ($\delta$)}"
		fname_stub = "avg_kend"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
			#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = "end of season" #L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot(0:(1/num_steps):1, kend[:,r,1], label=curr_label, color=col[r])
			end
			#legend(bbox_to_anchor=[.65,.95],loc="upper left", title=legendtitlestring) 
			legend(bbox_to_anchor=[0,.95],loc="upper left", title=legendtitlestring) 
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx)),
											yticks=(Array(miny:incy:maxy)),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
											#yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
											legend=:bottomright,
											#legend=:best,
											grid=false);
			for r = 1:num_rankings
				#cutoff_game = breakpoint_list[r]
				curr_label = ""
				if breakpoint_list[r] == 1	
					#curr_label = latexstring("$tmp", "\\mbox{ of season}")
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot!(0:(1/num_steps):1, kend[:,r,1], label=curr_label, linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot avg_games_tanked (# games tanked by draft ranking breakpoint)
		print("Plotting avg_games_tanked: average number of games tanked\n")
		miny = Int(floor(findmin(games_tanked[:,:,1])[1]))
		incy = 50 #Int(ceil((maxy - miny) / (5 * 10)) * 10)
		maxy = Int(ceil(findmax(games_tanked[:,:,1])[1] / incy) * incy)  #Int(floor(findmax(avg_games_tanked)[1]))
		titlestring = L"\mbox{Total games tanked}"
		xlabelstring = L"\mbox{Probability of tanking once eliminated}"
		ylabelstring = L"\mbox{Number of tanked games}"
		fname_stub = "avg_games_tanked"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
			#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot(0:(1/num_steps):1, games_tanked[:,r,1], label=curr_label, color=col[r])
			end
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx)),
											yticks=(Array(miny:incy:maxy)),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
											#yticks=(Array(miny:incy:maxy),[@sprintf("\$%.0f\$", i) for i in miny:incy:maxy]),
											legend=:topleft,
											#legend=:best,
											grid=false);
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot!(0:(1/num_steps):1, games_tanked[:,r,1], label=curr_label, linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot avg_already_tank
		print("Plotting avg_already_tank: average number of tanking teams\n")
		minx = 0
		incx = 0.1
		maxx = 1
		miny = Int(floor(findmin(already_tank[:,:,1])[1]))
		incy = 1
		maxy = Int(ceil(findmax(already_tank[:,:,1])[1]))
		titlestring = L"\mbox{Number of tanking teams}"
		xlabelstring = L"\mbox{Probability of tanking once eliminated}"
		ylabelstring = L"\mbox{Average number of tanking teams}"
		fname_stub = "avg_already_tank"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
			#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot(0:(1/num_steps):1, already_tank[:,r,1], label=curr_label, color=col[r])
			end
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx)),
											yticks=(Array(miny:incy:maxy)),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
											#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
											#yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
											legend=:topleft,
											grid=false);
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
				plot!(0:(1/num_steps):1, already_tank[:,r,1], label=curr_label, linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot avg_eliminated
		print("Plotting avg_eliminated: average number of eliminated teams by every game of the season\n")
		minx = 1
		maxx = num_games
		incx = (maxx - minx) / 5
		miny = 0 #Int(floor(findmin(math_eliminated[:,:,1])[1]))
		incy = 1
		maxy = num_teams - num_playoff_teams # Int(ceil(findmax(math_eliminated[:,:,1])[1]))
		titlestring = L"\mbox{Number of teams eliminated over time}"
		xlabelstring = L"\mbox{Percent of season elapsed}"
		ylabelstring = L"\mbox{Number of teams eliminated}"
		fname_stub = "avg_eliminated"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx), [@sprintf("%.0f", (100*i/num_games)) for i in minx:incx:maxx])
			yticks(Array(miny:incy:maxy))
			#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
			x = 1:num_games
			y = sum(math_eliminated[:,:,1], dims=1)[1,:] / (num_steps + 1)
			bar(x,y)
			#plot(x,y)
			#fill_between(x,y)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											xticks=(Array(minx:incx:maxx),[@sprintf("\$%.0f\$", (100*i/num_games)) for i in minx:incx:maxx]),
											yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
											legend=:none,
											grid=false);
			plot!(1:num_games, sum(math_eliminated[:,:,1], dims=1)[1,:] / (num_steps + 1), linetype=:bar);
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end
	end # if do_plotting
	return
end; # main_simulate

function main_parse(;do_plotting=true, mode=MODE, data_dir="../data", results_dir="../results")
	set_mode(mode)

	num_teams_eliminated_1314, num_games_tanked_1314, stats1314, critical_game1314 = parseNBASeason("games1314.csv", breakpoint_list, data_dir)
	num_teams_eliminated_1415, num_games_tanked_1415, stats1415, critical_game1415 = parseNBASeason("games1415.csv", breakpoint_list, data_dir)
	num_teams_eliminated_1516, num_games_tanked_1516, stats1516, critical_game1516 = parseNBASeason("games1516.csv", breakpoint_list, data_dir)
	num_teams_eliminated_1617, num_games_tanked_1617, stats1617, critical_game1617 = parseNBASeason("games1617.csv", breakpoint_list, data_dir)
	num_teams_eliminated_1718, num_games_tanked_1718, stats1718, critical_game1718 = parseNBASeason("games1718.csv", breakpoint_list, data_dir)

	# Retrieve data for avg_eliminated
	avg_eliminated = readdlm(string(results_dir, "/avg_eliminated", csvext), ',')
	num_steps = size(avg_eliminated)[1] - 1
	avg_eliminated = sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1)

	if false
		print("Year 2013-2014\n")
		print("num teams eliminated: ",num_teams_eliminated_1314,"\n")
		print("num possible games tanked: ",num_games_tanked_1314,"\n")

		print("\nYear 2014-2015\n")
		print("num teams eliminated: ",num_teams_eliminated_1415,"\n")
		print("num possible games tanked: ",num_games_tanked_1415,"\n")

		print("\nYear 2015-2016\n")
		print("num teams eliminated: ",num_teams_eliminated_1516,"\n")
		print("num possible games tanked: ",num_games_tanked_1516,"\n")

		print("\nYear 2016-2017\n")
		print("num teams eliminated: ",num_teams_eliminated_1617,"\n")
		print("num possible games tanked: ",num_games_tanked_1617,"\n")

		print("\nYear 2017-2018\n")
		print("num teams eliminated: ",num_teams_eliminated_1718,"\n")
		print("num possible games tanked: ",num_games_tanked_1718,"\n")
	end

	if (do_plotting)
		ind = [3,5,6] # needs to be ascending
		@assert ( length(breakpoint_list) in ind )
		num_years=5
		labels = [L"2013-2014", L"2014-2015", L"2015-2016", L"2016-2017", L"2017-2018"]
		col_labels = ["red", "orange", "green", "blue", "violet"]

		## Plot # games tanked
		print("Plotting num_games_tanked: number of games (possibly) tanked by the breakpoint mark\n")
		#num_games_tanked = zeros(Int, num_years, length(breakpoint_list))
		#num_games_tanked[1,:] = num_games_tanked_1314
		#num_games_tanked[2,:] = num_games_tanked_1415
		#num_games_tanked[3,:] = num_games_tanked_1516
		#num_games_tanked[4,:] = num_games_tanked_1617
		#num_games_tanked[5,:] = num_games_tanked_1718
		num_games_tanked = hcat(num_games_tanked_1314, num_games_tanked_1415, num_games_tanked_1516, num_games_tanked_1617, num_games_tanked_1718)
		num_games_tanked = num_games_tanked'
		#print(num_games_tanked,"\n")

		num_games_tanked_stacked = zeros(Int, num_years, length(ind))
		for tmp_i in 1:length(ind)
			i = ind[tmp_i]
			if tmp_i == 1
				num_games_tanked_stacked[:,tmp_i] = num_games_tanked[:,i]
			else
				prev_i = ind[tmp_i-1]
				num_games_tanked_stacked[:,tmp_i] = num_games_tanked[:,i] - num_games_tanked[:,prev_i]
			end
		end
		#print("num games tanked stacked: ",num_games_tanked_stacked,"\n")

		miny = 0 #Int(floor(findmin(num_games_tanked)[1]))
		maxy = Int(ceil(findmax(num_games_tanked)[1]))
		incy = 50 #(maxy - miny) / 5
		titlestring = L"\mbox{Number of games that could be tanked}"
		xlabelstring = L"\mbox{Season}"
		ylabelstring = L"\mbox{Number of possibly tanked games}"
		#legendtitlestring = L"\mbox{Draft ranking breakpoint}"
		legendtitlestring = L"\mbox{Breakpoint ($\delta$)}"
		fname_stub = "nba_num_teams_eliminated"
		fname_stub = "nba_num_games_tanked"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(1:num_years,["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]) 
			yticks(miny:incy:maxy)
			width = 0.75
			cumsum = zeros(Int, num_years, 1)
			for i in 1:length(ind)
				r = ind[i]
				curr_label = ""
				if breakpoint_list[r] == 1
					curr_label=L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end	
				if i == 1
					bar(1:num_years, num_games_tanked_stacked[:,i], label=curr_label, width = width)
				else
					bar(1:num_years, num_games_tanked_stacked[:,i], bottom=cumsum[:,1], label=curr_label, width = width)
				end
				cumsum += num_games_tanked_stacked[:,i]
			end
			#legend(loc="best", title=legendtitlestring)
			legend(bbox_to_anchor=[1,.9],loc="upper right", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			ctg = repeat(vcat([latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}") for r in ind if r < length(breakpoint_list)], L"\mbox{end of season}"), inner=num_years)
			fig = groupedbar(num_games_tanked_stacked,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
											yticks=(Array(miny:incy:maxy),[@sprintf("\$%d\$", i) for i in miny:incy:maxy]),
											group=ctg,
											lw=0,
											bar_position=:stack,
											legend=:best,
											grid=false,
											show=false)
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot number of teams eliminated
		print("Plotting num_teams_eliminated: number of eliminated teams by every game of the season\n")
		num_games = 15 * 82;
		num_teams_eliminated = hcat(num_teams_eliminated_1314, num_teams_eliminated_1415, num_teams_eliminated_1516, num_teams_eliminated_1617, num_teams_eliminated_1718)
		num_teams_eliminated = num_teams_eliminated'

		minx = 1 / num_games
		maxx = num_games / num_games
		incx = (maxx - minx) / 5
		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
		incy = 1
		titlestring = L"\mbox{Number of teams eliminated over time}"
		xlabelstring = L"\mbox{Percent of season elapsed}"
		ylabelstring = L"\mbox{Number of teams eliminated}"
		legendtitlestring = L"\mbox{Season}"
		fname_stub = "nba_num_teams_eliminated"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx), [@sprintf("%.0f", (100*i/maxx)) for i in minx:incx:maxx])
			#xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
			for s = 1:size(num_teams_eliminated)[1]
				curr_num_games = length(num_teams_eliminated[s,:])
				curr_label = labels[s]
				plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), num_teams_eliminated[s,:], label=curr_label, color=col_labels[s]);
			end
			curr_num_games = length(avg_eliminated)
			plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), avg_eliminated, label="simulated", color="black", linestyle="dashed")
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx),[@sprintf("\$%.0f\$", (100*i/num_games)) for i in minx:incx:maxx]),
											yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
											legend=:topleft,
											grid=false);
			for s = 1:size(num_teams_eliminated)[1]
				curr_num_games = length(num_teams_eliminated[s,:])
				curr_label = labels[s]
				plot!(1:curr_num_games, num_teams_eliminated[s,:], label=curr_label, linecolor=col_labels[s]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		# Print number of teams eliminated
#		print("Plotting num_teams_eliminated: number of eliminated teams by the breakpoint mark\n")
#		#num_teams_eliminated = zeros(Int, num_years, length(breakpoint_list))
#		#num_teams_eliminated[1,:] = num_teams_eliminated_1314
#		#num_teams_eliminated[2,:] = num_teams_eliminated_1415
#		#num_teams_eliminated[3,:] = num_teams_eliminated_1516
#		#num_teams_eliminated[4,:] = num_teams_eliminated_1617
#		#num_teams_eliminated[5,:] = num_teams_eliminated_1718
#		num_teams_eliminated = hcat(num_teams_eliminated_1314, num_teams_eliminated_1415, num_teams_eliminated_1516, num_teams_eliminated_1617, num_teams_eliminated_1718)
#		num_teams_eliminated = num_teams_eliminated'
#		#print(num_teams_eliminated,"\n")
#
#		num_teams_eliminated_stacked = zeros(Int, num_years, length(breakpoint_list))
#		num_teams_eliminated_stacked[:,1] = num_teams_eliminated[:,1]
#		num_teams_eliminated_stacked[:,2] = num_teams_eliminated[:,2] - num_teams_eliminated[:,1]
#		num_teams_eliminated_stacked[:,3] = num_teams_eliminated[:,3] - num_teams_eliminated[:,2]
#		num_teams_eliminated_stacked[:,4] = num_teams_eliminated[:,4] - num_teams_eliminated[:,3]
#		num_teams_eliminated_stacked[:,5] = num_teams_eliminated[:,5] - num_teams_eliminated[:,4]
#
#		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
#		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
#		inc = floor((maxy - miny) / 5)
#		ctg = repeat(vcat([latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}") for r in ind if r < length(breakpoint_list)], L"\mbox{end of season}"), inner=num_years)
#		fig = groupedbar(num_teams_eliminated[:,ind], 
#										xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
#										yticks=(Array(miny:inc:maxy),[@sprintf("\$%d\$", i) for i in miny:inc:maxy]),
#										group=ctg,
#										lw=0,
#										bar_position=:dodge,
#										#bar_position=:stack,
#										title=L"\mbox{Number of teams eliminated by breakpoint mark}",
#										xlab=L"\mbox{Season}", 
#										ylab=L"\mbox{Number of teams eliminated}",
#										legend=:best,
#										legendfont=6,
#										legendtitle=L"\mbox{Draft ranking breakpoint}",
#										titlefont=12,
#										tickfont=8,
#										grid=false,
#										show=false)
#		savefig(fig, string(results_dir,"/nba_num_teams_eliminated",ext));
	end

	return
end; # main_parse

function rankings_are_noisy(;do_simulation=true, num_replications=1000, do_plotting=true, mode=MODE, results_dir="../results")
	set_mode(mode)

	## Variables that need to be set
	num_rounds_set = [1 2 3 4 5 10 100 1000];
	num_steps = 50; # number of subdivisions of [0.5,1]
	prob = 0.5:0.5/num_steps:1;

	## Set constants
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2);

	## For output
	avg_kend = zeros(Float64, num_steps+1, length(num_rounds_set));

	if do_simulation
		for step_ind = 1:num_steps+1
			print("Step ", step_ind, " of ", num_steps+1, "\n")
			gamma = prob[step_ind];
			for num_rounds_ind in 1:length(num_rounds_set)
				num_rounds = num_rounds_set[num_rounds_ind];
				for rep = 1:num_replications
					stats = zeros(Int, num_teams, 2);
					for i = 1:num_teams
						stats[i,1] = i;
						stats[i,2] = 0;
					end

					for round_ind = 1:num_rounds
						for i = 1:num_teams
							for j = i+1:num_teams
								team_i_wins = teamWillWinNoTanking(i,j,gamma,true_strength,mode)

								# Do updates
								for k in [i,j]
									team_k_wins = (k == i) ? team_i_wins : !team_i_wins
									stats[k,2] = stats[k,2] + team_k_wins
								end
							end
						end
					end # loop over rounds

					## Calculate Kendell tau distance for this round
					avg_kend[step_ind, num_rounds_ind] = avg_kend[step_ind, num_rounds_ind] + kendtau(stats, 2, true_strength, mode, num_teams - num_playoff_teams) / num_replications
				end # loop over replications
			end # loop over num_rounds_set
		end # loop over steps
		writedlm(string(results_dir, "/noisy_ranking", csvext), avg_kend, ',')
	else
		avg_kend = readdlm(string(results_dir, "/noisy_ranking", csvext), ',')
	end # if do_simulation

	if do_plotting
		## Plot noisy ranking
		print("Plotting noisy_ranking: shows low fidelity of ranking unless teams play each other many times\n")
		minx = 0.5
		incx = 0.1
		maxx = 1
		miny = Int(ceil(findmin(avg_kend)[1]));
		incy = 50 #Int(floor((maxy - miny) / 5))
		maxy = Int(floor(findmax(avg_kend)[1]));
		#maxy = Int(ceil(maxy/incy)*incy)
		titlestring = L"\mbox{Effect of number of rounds on ranking accuracy}"
		xlabelstring = L"\mbox{Probability better team wins}"
		ylabelstring = L"\mbox{Distance from true ranking of all teams}"
		legendtitlestring = L"\mbox{Rounds}"
		fname_stub = "noisy_ranking"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
			for r = 1:length(num_rounds_set)
				curr_label = latexstring(num_rounds_set[r])
				plot(prob, avg_kend[:,r], label=curr_label);
			end
			legend(loc="upper right", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)

			#size = fig[:get_size_inches]()
			#print("Which should result in a ", DPI*size[1], " x ", DPI*size[2], " image")

			close()
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx)),
											yticks=(Array(miny:incy:maxy)),
											#xticks=(Array(minx:incx:maxx),["\$$i\$" for i in minx:incx:maxx]),
											#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
											#legendtitle=L"\mbox{Rounds}",
											legend=:topright)
			for r = 1:length(num_rounds_set)
				curr_label = latexstring(num_rounds_set[r])
				plot!(prob, avg_kend[:,r], label=curr_label);
								#linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end
	end # if do_plotting
end; # rankings_are_noisy

function tanking_unit_tests()
	test_ranking = 1:num_teams
	test_strength = 1:num_teams
	@assert ( kendtau_sorted(test_ranking, test_strength, 1) == Int(num_teams * (num_teams-1) / 2) )
	test_strength = num_teams:-1:1
	@assert ( kendtau_sorted(test_ranking, test_strength, 1) == 0 )
end # tanking_unit_tests
