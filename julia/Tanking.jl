######################################
# On Tanking and Competitive Balance #
######################################
# Aleksandr M. Kazachkov
# Shai Vardi
###
# February 2020
###
# Implemented are the bilevel, Gold, Lenten, and NBA lottery-based rankings (both pre- and post-2019 changes)
# 1. Bilevel: set a breakpoint; rank nonplayoff teams based on their relative order at the breakpoint
# 2. Gold: rank based on number of wins after elimination
# 3. Lenten: reverse of order in which teams were eliminated
# 4. NBA lottery-based system: see odds below, which are applied to the end-of-season standings
###

module Tanking

## Required dependencies
import Random
using DelimitedFiles
using LaTeXStrings
using Printf
import Distributions.Beta
include("mathelim.jl") # needed for simulate.jl and utility.jl
include("utility.jl") # imports MODE definitions, needed for parse.jl and simulate.jl
include("BT.jl")
include("parse.jl")
include("simulate.jl")

using Combinatorics # for permutations

## NBA draft odds (for non-playoff teams, in reverse order)
# If teams are tied, then those teams receive odds that are the average of the odds for the positions they occupy
# Keep these in order that teams are ranked (not for the draft, but by wins)
nba_odds_old = [.250, .199, .156, .119, .088, .063, .043, .028, .017, .011, .008, .007, .006, .005]
nba_odds_new = [.140, .140, .140, .125, .105, .090, .075, .060, .045, .030, .020, .015, .010, .005]
nba_odds_flat = [1. / 14. for i in 1:14]
nba_odds_old = nba_odds_old[14:-1:1]
nba_odds_new = nba_odds_new[14:-1:1]
nba_odds_list = [nba_odds_new, nba_odds_old, nba_odds_flat]
nba_num_lottery = [4, 3, 4]
@assert( length(nba_num_lottery) == length(nba_odds_list) )

## When the (draft) ranking will be set as fraction of number games
#breakpoint_list = [4//8; 5//8; 6//8; 7//8; 1]
breakpoint_list = [1//2; 2//3; 3//4; 5//6; 7//8; 1]
num_rankings = length(breakpoint_list)
#shape = [:vline, :utriangle, :rect, :x, :triangle, :circle]
col   = ["red",   "orange", "green",   "blue",   "violet", "black", "gray"]
style = ["solid", "dashed", "dashdot", "dotted", "dashed", "solid", "solid"]
shape = ["+",  "s",   ".",    "",   "",   "",  ""]
shapesize = [5, 2, 5, 5, 5, 5, 5]
#color_for_cutoff_point = ["c" "b" "m" "r" "k"]
num_teams = 30 # number of teams
num_playoff_teams = Int(2^ceil(log(2, num_teams / 2)))

## Set ranking type
# 1: strict: teams are strictly ordered 1 \succ 2 \succ \cdots \succ 30
# 2: ties: [1,5] \succ [6,10] \succ \cdots \succ [26,30]
# 3: BT_distribution: Bradley-Terry with P(i>j) = p_i / (p_i + p_j); must also set distribution, where default is each team gets a strength score from U[0,1]
# 4: BT_exponential: Bradley-Terry with P(i>j) = exp(p_i) / (exp(p_i) + exp(p_j)); must also set distribution, where default is each team gets a strength score from U[0,1]; can consider others such as, e.g., using Beta(alpha=3, beta=5)
MODE = STRICT # imported from utility.jl

ranking_type = ""
true_strength = []
csvext = ".csv"
ext_folder = "pdf"
lowext_folder = "png"
ext = string(ranking_type,".",ext_folder)
lowext = string(ranking_type,"_low",".",lowext_folder)

"""
    set_mode

Set global variables based on which mode is being used
"""
function set_mode(mode=MODE, extra=nothing)
  tmp_true_strength = 0
	if mode == NONE
		cssvext = ".csv"
	elseif mode == STRICT
		# 1 \succ 2 \succ \cdots \succ 30
		global ranking_type="_strict"
		tmp_true_strength = num_teams:-1:1
	elseif mode == TIES
		# [1,5] \succ [6,10] \succ \cdots \succ [26,30]
		global ranking_type="_ties"
		tmp_true_strength = [Int(ceil(i/5)) for i in num_teams:-1:1] # allows for ties # old: [i:i+4 for i in 1:5:num_teams-4]
	elseif mode == BT_DISTR
		# Options to consider:
		# uniform distribution (same as beta(1,1))
		# --> essentially perfectly imbalanced
		# nonuniform distribution, beta(2,2), or maybe we should do some kind of bimodal distribution
		# --> more weight on middle teams, smaller probability of very weak or very strong teams
		global ranking_type="_BT_distr"
		#tmp_true_strength = rand(30)
    distr = Beta(2,5)
    num_repeats = 1000
    tmp_true_strength = sort(rand(distr, 30), rev=true)
    for i = 2:num_repeats
      tmp_true_strength += sort(rand(distr, 30), rev=true)
    end
    tmp_true_strength /= num_repeats
	elseif mode == BT_EXPONENTIAL
		global ranking_type="_BT_exp"
		#tmp_true_strength = rand(30)
    distr = Beta(2,5)
    tmp_true_strength = rand(distr, 30)
  elseif mode == BT_ESTIMATED
    global ranking_type="_BT_est"
    # Need to estimate true_strength
    #BT_avg() BT_MLE()
    tmp_true_strength = BT_MLE()
	end

  ## Append extra string to ranking_type, if provided
  if isnothing(extra)
    # Do nothing
  elseif isa(extra,String)
    extra = filter(x -> !isspace(x), extra) # remove spaces
    global ranking_type = string(ranking_type,extra)
  elseif isa(extra,UnitRange)
    global ranking_type = string(ranking_type,extra)
  elseif isa(extra,Array) && length(extra) > 0
    for i = 1:length(extra)
      if !isa(extra[i], Number)
        error("Entries of `extra` array need to be numbers, but value of array at index $i = ", extra[i])
      end
    end
    # Check if we can make it a unit range
    is_range = length(extra) > 1 ? true : false
    for i = 2:length(extra)
      if extra[i] != extra[i-1]+1
        is_range = false
        break
      end
    end
    if is_range
      extra = UnitRange(extra[1],extra[length(extra)])
    end
    extra = filter(x -> !isspace(x), string(extra)) # remove spaces
    global ranking_type = string(ranking_type,extra)
  end # append extra string

	global csvext = string(ranking_type,".csv")
	global ext = string(ranking_type,".",ext_folder)
	global lowext = string(ranking_type,"_low",".",lowext_folder)
  global true_strength = sort(tmp_true_strength, rev=true)
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
    pygui(:default)
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

"""
    clean_selected_steps

Ensure selected steps is in the proper format (an `Int`, `Int` array, or a `UnitRange`) 
"""
function clean_selected_steps(selected_steps = nothing)
  if isnothing(selected_steps) || isa(selected_steps, Int) || isa(selected_steps, UnitRange)
    return selected_steps
  elseif !isa(selected_steps, Array)
    error("selected_steps needs to be an Int, Int array, or UnitRange; provided is ", selected_steps)
  else
    # We have an array, but we should make sure each entry is an Int
    all_ok = true
    for i = 1:length(selected_steps)
      if !isa(selected_steps[i], Int) && !isa(selected_steps[i], UnitRange)
        error("selected_steps needs to be an Int, Int array, or UnitRange, but entry $i is ", selected_steps[i])
      elseif !isa(selected_steps[i], Int)
        all_ok = false
      end
    end
    if all_ok
      return selected_steps
    end

    # Otherwise we need to do some fixing
    x = Array{Any}(selected_steps)
    i = 1; len = length(x)
    while i <= len
      curr_len = length(x[i])
      splice!(x, i, x[i])
      i += curr_len
      len += curr_len - 1
    end
    return x
  end
end # clean_selected_steps

"""
    main_simulate

Simulate a season and plot output

Parameters
---
  * `do_simulation`: when false, read data from files in results_dir
  * `num_replications`: how many times to simulate each data point
  * `do_plotting`: if false, only gather data, without plotting it
  * `mode`: which kind of true ranking is used;
       1 or 2: 
         When two non-tanking teams or two tanking teams play each other, 
         the better team wins with probability gamma.
         When a tanking team plays a non-tanking team, the tanking team always loses.
       3 or 4:
         Variants of (Zermelo-)Bradley-Terry model used to determine who wins each game.

  * `results_dir`: where results should be saved and can be found
  * `num_rounds`: a round consists of each team playing each other team
  * `num_steps`: discretization of [0,1] for tanking probability
  * `gamma`: probability a better-ranked team wins over a worse-ranked team
  * `math_elim_mode`: identify when a team is eliminated, and will stop tanking;
    0: use effective elimination, 
    1: use mathematical elimination, but calculated by heuristics only, 
    2: use math elim, binary MIP, team-wise formulation,
    3: use math elim, general integer MIP, team-wise formulation,
    4: use math elim, binary MIP, cutoff formulation,
    5: use math elim, general integer MIP, cutoff formulation,
    <0: use effective elimination for tanking, but calculate mathematical elimination.
  * `selected_steps`: subset of steps we wish to actually get results for; this should be an `Int`, `Int` array, or a `UnitRange`
"""
function main_simulate(;do_simulation = true, num_replications = 100000, 
    do_plotting = true, mode = MODE, results_dir = "../results", 
    num_rounds = 3, num_steps = num_teams, gamma = 0.71425, 
    math_elim_mode = -2, selected_steps = nothing)
  Random.seed!(628) # for reproducibility
  selected_steps = clean_selected_steps(selected_steps)
	set_mode(mode, selected_steps)

	## Variables that need to be set
	## end variables that need to be set

	## Set constants
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
	num_games           = num_rounds * num_games_per_round

	## For output
	kend                  = 0 # [step,cutoff,stat], holds KT distance for each tanking probability and cutoff for draft ranking
  kend_nba              = 0 # [step,odds,stat], KT distance based on NBA ranking
	kend_gold             = 0 # [:,stat], ranking based on number of wins since elimination point
	kend_lenten           = 0 # [step,stat], ranking in order of first-to-eliminated
	games_tanked          = 0 # [step,cutoff,stat], number games tanked by cutoff
	already_tank          = 0 # [step,cutoff,stat], number teams already tanking by the cutoff
	math_eliminated       = 0 # [step,game,stat], number teams eliminated by each game
	eff_eliminated        = 0 # [step,game,stat], number teams eliminated by each game
  num_mips              = 0 # [step,stat], number of MIPs solved
  num_unelim            = 0 # [step,stat], number of teams eff elim, then uneliminated
  avg_rank_strat        = 0 # [step,stat], average rank of strategic teams
  avg_rank_moral        = 0 # [step,stat], average rank of moral teams
  avg_elim_rank_strat   = 0 # [step,stat], average rank of eliminated strategic teams
  avg_elim_rank_moral   = 0 # [step,stat], average rank of eliminated moral teams
  avg_diff_rank_strat   = 0 # [step,stat], average of differences between true and calculated ranks of strategic teams
  avg_diff_rank_moral   = 0 # [step,stat], average of differences between true and calculated ranks of moral teams
  num_missing_case      = 0 # [step,stat], number of times incomplete case is encountered (\delta < \tau_i \le t + conditions)

  ## Stats we keep
  avg_stat    = 1
  stddev_stat = 2
  min_stat    = 3
  max_stat    = 4
  num_stats   = 4
  prefix      = ["avg_", "stddev_", "min_", "max_"]

	## Do simulation or retrieve data
	if do_simulation
		## Do simulation
    kend, kend_nba, kend_gold, kend_lenten, games_tanked, already_tank, 
      math_eliminated, eff_eliminated, num_mips, num_unelim, 
      avg_rank_strat, avg_rank_moral,
      avg_elim_rank_strat, avg_elim_rank_moral,
      avg_diff_rank_strat, avg_diff_rank_moral,
      num_missing_case = 
        simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, gamma, breakpoint_list, nba_odds_list, nba_num_lottery, true_strength, mode, math_elim_mode, selected_steps, false)

    for stat = 1:num_stats
      writedlm(string(results_dir, "/", prefix[stat], "kend", csvext), kend[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "kend_lenten", csvext), kend_lenten[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "kend_nba", csvext), kend_nba[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "games_tanked", csvext), games_tanked[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "already_tank", csvext), already_tank[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "math_eliminated", csvext), math_eliminated[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "eff_eliminated", csvext), eff_eliminated[:,:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "num_mips", csvext), num_mips[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "num_unelim", csvext), num_unelim[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_rank_strat", csvext), avg_rank_strat[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_rank_moral", csvext), avg_rank_moral[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_elim_rank_strat", csvext), avg_elim_rank_strat[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_elim_rank_moral", csvext), avg_elim_rank_moral[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_diff_rank_strat", csvext), avg_diff_rank_strat[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "avg_diff_rank_moral", csvext), avg_diff_rank_moral[:,stat], ',')
      writedlm(string(results_dir, "/", prefix[stat], "num_missing_case", csvext), num_missing_case[:,stat], ',')
    end
    writedlm(string(results_dir, "/", "kend_gold", csvext), kend_gold, ',')
	else
    ## Resize things
    num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
    num_games_total     = num_rounds * num_games_per_round
    kend                = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    kend_nba            = zeros(Float64, num_steps+1, length(nba_odds_list), num_stats)
    kend_gold           = zeros(Float64, 1, num_stats)
    kend_lenten         = zeros(Float64, num_steps+1, num_stats)
    games_tanked        = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    already_tank        = zeros(Float64, num_steps+1, length(breakpoint_list), num_stats)
    math_eliminated     = zeros(Float64, num_steps+1, num_games_total, num_stats)
    eff_eliminated      = zeros(Float64, num_steps+1, num_games_total, num_stats)
    num_mips            = zeros(Float64, num_steps+1, num_stats)
    num_unelim          = zeros(Float64, num_steps+1, num_stats)
    avg_rank_strat      = zeros(Float64, num_steps+1, num_stats)
    avg_rank_moral      = zeros(Float64, num_steps+1, num_stats)
    avg_elim_rank_strat = zeros(Float64, num_steps+1, num_stats)
    avg_elim_rank_moral = zeros(Float64, num_steps+1, num_stats)
    avg_diff_rank_strat = zeros(Float64, num_steps+1, num_stats)
    avg_diff_rank_moral = zeros(Float64, num_steps+1, num_stats)
    num_missing_case    = zeros(Float64, num_steps+1, num_stats)
    for stat = 1:num_stats
      kend[:,:,stat]            = readdlm(string(results_dir, "/", prefix[stat], "kend", csvext), ',')
      kend_lenten[:,stat]       = readdlm(string(results_dir, "/", prefix[stat], "kend_lenten", csvext), ',')
      kend_nba[:,:,stat]        = readdlm(string(results_dir, "/", prefix[stat], "kend_nba", csvext), ',')
      games_tanked[:,:,stat]    = readdlm(string(results_dir, "/", prefix[stat], "games_tanked", csvext), ',')
      already_tank[:,:,stat]    = readdlm(string(results_dir, "/", prefix[stat], "already_tank", csvext), ',')
      math_eliminated[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], "math_eliminated", csvext), ',')
      eff_eliminated[:,:,stat]  = readdlm(string(results_dir, "/", prefix[stat], "eff_eliminated", csvext), ',')
      num_mips[:,stat]          = readdlm(string(results_dir, "/", prefix[stat], "num_mips", csvext), ',')
      num_unelim[:,stat]        = readdlm(string(results_dir, "/", prefix[stat], "num_unelim", csvext), ',')
      avg_rank_strat[:,stat]    = readdlm(string(results_dir, "/", prefix[stat], "avg_rank_strat", csvext), ',') 
      avg_rank_moral[:,stat]    = readdlm(string(results_dir, "/", prefix[stat], "avg_rank_moral", csvext), ',') 
      avg_elim_rank_strat[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "avg_elim_rank_strat", csvext), ',') 
      avg_elim_rank_moral[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "avg_elim_rank_moral", csvext), ',') 
      avg_diff_rank_strat[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "avg_diff_rank_strat", csvext), ',') 
      avg_diff_rank_moral[:,stat] = readdlm(string(results_dir, "/", prefix[stat], "avg_diff_rank_moral", csvext), ',') 
      num_missing_case[:,stat]  = readdlm(string(results_dir, "/", prefix[stat], "num_missing_case", csvext), ',')
    end
    kend_gold = readdlm(string(results_dir, "/", "kend_gold", csvext), ',')
	end
  num_eliminated = (math_elim_mode > 0) ? math_eliminated : eff_eliminated

	if (do_plotting)
		## Plot avg_kend (Kendell tau distance) for the bilevel ranking
		print("Plotting avg_kend: average swap distance\n")
    if num_steps == num_teams
      minx = 0
      incx = 5
      maxx = num_teams
    else
      minx = 0
      incx = 0.1
      maxx = 1
    end
		miny = Int(floor(findmin(kend[:,:,1])[1]))
		incy = 1
		maxy = 31 #Int(ceil(findmax(kend[:,:,1])[1]))
		titlestring = L"\mbox{Effect of $\delta$ on bilevel ranking of non-playoff teams}"
    if num_steps == num_teams
      xlabelstring = L"\mbox{Number of selfish teams}"
    else
      xlabelstring = L"\mbox{Probability of tanking once eliminated}"
    end
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
        if num_steps == num_teams
          plot(0:1:num_teams, kend[:,r,avg_stat], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        else
				  plot(0:(1/num_steps):1, kend[:,r,avg_stat], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        end
			end

      # Also plot kend_nba (old and new) and kend_lenten (flat line)
      curr_label = "NBA (old)"
      if num_steps == num_teams
        plot(0:1:num_teams, kend_nba[:,1,avg_stat], label=curr_label, color="gray")
      else
        plot(0:(1/num_steps):1, kend_nba[:,1,avg_stat], label=curr_label, color="gray")
      end
      curr_label = "NBA (new)"
      if num_steps == num_teams
        plot(0:1:num_teams, kend_nba[:,2,avg_stat], label=curr_label, color="gray", marker=".", markersize=5)
      else
        plot(0:(1/num_steps):1, kend_nba[:,2,avg_stat], label=curr_label, color="gray", marker=".", markersize=5)
      end
      curr_label = "Lenten"
      if num_steps == num_teams
        plot(0:1:num_teams, kend_lenten[:,avg_stat], label=curr_label, color="gray", marker="x", markersize=5)
      else
        plot(0:(1/num_steps):1, kend_lenten[:,avg_stat], label=curr_label, color="gray", marker="x", markersize=5)
      end
      #plot(0:(1/num_steps):1, [kend_lenten[avg_stat] for i in 0:(1/num_steps):1], label=curr_label, color="gray", marker="x")
      #axhline(kend_lenten[avg_stat], label=curr_label, color="gray", marker="x")

			#legend(bbox_to_anchor=[.65,.95],loc="upper left", title=legendtitlestring) 
			#legend(bbox_to_anchor=[0,.95],loc="upper left", title=legendtitlestring) 
			legend(bbox_to_anchor=[1.0, 0.5],loc="center left", title=legendtitlestring) 
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
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
		miny = 0 #Int(floor(findmin(games_tanked[:,:,1])[1]))
		incy = 50 #Int(ceil((maxy - miny) / (5 * 10)) * 10)
		maxy = Int(ceil(findmax(games_tanked[:,:,1])[1] / incy) * incy)  #Int(floor(findmax(avg_games_tanked)[1]))
		titlestring = L"\mbox{Total games tanked}"
    if num_steps == num_teams
      xlabelstring = L"\mbox{Number of selfish teams}"
    else
      xlabelstring = L"\mbox{Probability of tanking once eliminated}"
    end
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
      ylim(ymin = -5)
      ylim(ymax = maxy)
			#yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
			for r = 1:num_rankings
				curr_label = ""
				if breakpoint_list[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}")
				end
        if num_steps == num_teams
          plot(0:1:num_teams, games_tanked[:,r,1], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        else
          plot(0:(1/num_steps):1, games_tanked[:,r,1], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        end
			end
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
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
    if num_steps == num_teams
      minx = 0
      incx = 5
      maxx = num_teams
    else
      minx = 0
      incx = 0.1
      maxx = 1
    end
    miny = Int(floor(findmin(already_tank[:,:,1])[1]))
		incy = 1
		maxy = Int(ceil(findmax(already_tank[:,:,1])[1]))
		titlestring = L"\mbox{Number of tanking teams}"
    if num_steps == num_teams
      xlabelstring = L"\mbox{Number of selfish teams}"
    else
      xlabelstring = L"\mbox{Probability of tanking once eliminated}"
    end
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
        if num_steps == num_teams
          plot(0:1:num_teams, already_tank[:,r,1], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        else
          plot(0:(1/num_steps):1, already_tank[:,r,1], label=curr_label, color=col[r], linestyle=style[r], marker=shape[r], markersize=shapesize[r])
        end
			end
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
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
		miny = 0 #Int(floor(findmin(num_eliminated[:,:,1])[1]))
		incy = 1
		maxy = num_teams - num_playoff_teams # Int(ceil(findmax(num_eliminated[:,:,1])[1]))
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
			y = sum(num_eliminated[:,:,1], dims=1)[1,:] / (num_steps + 1)
			bar(x,y)
			#plot(x,y)
			#fill_between(x,y)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											xticks=(Array(minx:incx:maxx),[@sprintf("\$%.0f\$", (100*i/num_games)) for i in minx:incx:maxx]),
											yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
											legend=:none,
											grid=false);
			plot!(1:num_games, sum(num_eliminated[:,:,1], dims=1)[1,:] / (num_steps + 1), linetype=:bar);
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

    ## Plot avg_rank (strat vs moral)
		print("Plotting avg_rank: average rank of strategic vs moral teams\n")
    if num_steps == num_teams
      minx = 0
      incx = 5
      maxx = num_teams
    else
      minx = 0
      incx = 0.1
      maxx = 1
    end
    miny = Int(floor(findmin(avg_rank_moral[1:num_steps,1])[1]))
		maxy = Int(ceil(findmax(avg_rank_strat[:,1])[1]))
		incy = 1
		titlestring = L"\mbox{Average rank of selfish versus moral teams}"
    if num_steps == num_teams
      xlabelstring = L"\mbox{Number of selfish teams}"
    else
      xlabelstring = L"\mbox{Probability of tanking once eliminated}"
    end
		ylabelstring = L"\mbox{Average rank}"
		legendtitlestring = L"\mbox{Team Type}"
		fname_stub = "avg_rank"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
      if num_steps == num_teams
        curr_label = "Moral"
        plot(0:1:num_teams-1, avg_rank_moral[1:num_steps,avg_stat], label=curr_label, color="green", linestyle="solid")
        curr_label = "Selfish"
        plot(1:1:num_teams, avg_rank_strat[2:num_steps+1,avg_stat], label=curr_label, color="red", linestyle="dashed")
      else
        curr_label = "Moral"
        plot(0:(1/num_steps):1-1/num_steps, avg_rank_moral[1:num_steps,avg_stat], label=curr_label, color="green", linestyle="solid")
        curr_label = "Selfish"
        plot(1/num_steps:(1/num_steps):1, avg_rank_strat[2:num_steps+1,avg_stat], label=curr_label, color="red", linestyle="dashed")
      end
			legend(loc="upper right", title=legendtitlestring) 
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
    end

    ## Plot avg_elim_rank (strat vs moral)
		print("Plotting avg_elim_rank: average rank of eliminated strategic vs moral teams\n")
    if num_steps == num_teams
      minx = 0
      incx = 5
      maxx = num_teams
    else
      minx = 0
      incx = 0.1
      maxx = 1
    end
		miny = Int(floor(findmin(avg_elim_rank_moral[2:num_steps,1])[1]))
		maxy = Int(ceil(findmax(avg_elim_rank_strat[2:num_steps,1])[1]))
		incy = 1
		titlestring = L"\mbox{Average rank of selfish versus moral non-playoff teams}"
    if num_steps == num_teams
      xlabelstring = L"\mbox{Number of selfish teams}"
    else
      xlabelstring = L"\mbox{Probability of tanking once eliminated}"
    end
		ylabelstring = L"\mbox{Average rank}"
		legendtitlestring = L"\mbox{Team Type}"
		fname_stub = "avg_elim_rank"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

		if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
			xlabel(xlabelstring)
			ylabel(ylabelstring)
			xticks(Array(minx:incx:maxx))
			yticks(Array(miny:incy:maxy))
      if num_steps == num_teams
        curr_label = "Moral"
        plot(1:1:num_teams-1, avg_elim_rank_moral[2:num_steps,avg_stat], label=curr_label, color="green", linestyle="solid")
        curr_label = "Selfish"
        plot(1:1:num_teams-1, avg_elim_rank_strat[2:num_steps,avg_stat], label=curr_label, color="red", linestyle="dashed")
      else
        curr_label = "Moral"
        plot(1/num_steps:(1/num_steps):1-1/num_steps, avg_elim_rank_moral[2:num_steps,avg_stat], label=curr_label, color="green", linestyle="solid")
        curr_label = "Selfish"
        plot(1/num_steps:(1/num_steps):1-1/num_steps, avg_elim_rank_strat[2:num_steps,avg_stat], label=curr_label, color="red", linestyle="dashed")
      end
			legend(loc="upper right", title=legendtitlestring) 
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
    end
	end # if do_plotting
	return
end; # main_simulate

"""
    main_parse

Parse data from 2004-2019, except 2011-12 (lockout year)
"""
function main_parse(;do_plotting=true, mode=MODE, data_dir="../data", results_dir="../results")
  Random.seed!(628) # for reproducibility
	set_mode(mode)

  #years = ["games1314.csv", "games1415.csv", "games1516.csv", "games1617.csv", "games1718.csv", "games1819.csv"]
  years = ["games0405.csv", "games0506.csv", "games0607.csv", "games0708.csv", "games0809.csv", "games0910.csv", "games1011.csv", "games1213.csv", "games1314.csv", "games1415.csv", "games1516.csv", "games1617.csv", "games1718.csv", "games1819.csv"]
  num_years = length(years)

  num_teams = 30
  num_games_per_team = 82
  num_games_total = Int(num_games_per_team * num_teams / 2)
  size_of_stats = 6
  num_teams_eliminated = zeros(Float64, num_years, num_games_total)
  num_games_tanked = zeros(Int, num_years, length(breakpoint_list))
  #stats = zeros(Float64, num_years, num_teams, size_of_stats)
  #critical_game = zeros(Int, num_years, num_teams, 3)

  for yr = 1:num_years
    num_teams_eliminated[yr,:], num_games_tanked[yr,:], _, _, _ = parseNBASeason(years[yr], breakpoint_list, data_dir)
  end # loop over years

	# Retrieve data for avg_eliminated
	avg_eliminated = readdlm(string(results_dir, "/avg_eff_eliminated", csvext), ',')
	num_steps = size(avg_eliminated)[1] - 1
	avg_eliminated = sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1)

	if (do_plotting)
		ind = [3,5,6] # needs to be ascending
		@assert ( length(breakpoint_list) in ind )
		#labels = [L"2013-2014", L"2014-2015", L"2015-2016", L"2016-2017", L"2017-2018"]
		#col_labels = ["red", "orange", "green", "blue", "violet"]
    #labels = [L"2004-05", L"2005-06", L"2006-07", L"2007-08", L"2008-09", L"2009-10", L"2010-11", L"2012-13", L"2013-14", L"2014-15", L"2015-16", L"2016-17", L"2017-18", L"2018-19"]
    col_labels = []
    labels = [L"04-05", L"05-06", L"06-07", L"07-08", L"08-09", L"09-10", L"10-11", L"12-13", L"13-14", L"14-15", L"15-16", L"16-17", L"17-18", L"18-19"]
    @assert ( length(labels) == num_years )
    @assert ( (length(col_labels) == 0) || (length(col_labels) == num_years) )

		## Plot # games tanked
		print("Plotting num_games_tanked: number of games (possibly) tanked by the breakpoint mark\n")
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

		miny = 0 #Int(floor(findmin(num_games_tanked)[1]))
		incy = 50 #(maxy - miny) / 5
		maxy = incy * Int(ceil(findmax(num_games_tanked)[1] / incy))
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
			#xticks(1:num_years,["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]) 
			xticks(1:num_years,labels)
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
			close(fig)
		else
			ctg = repeat(vcat([latexstring(numerator(breakpoint_list[r]),"/",denominator(breakpoint_list[r]), "\\mbox{ of season}") for r in ind if r < length(breakpoint_list)], L"\mbox{end of season}"), inner=num_years)
			fig = groupedbar(num_games_tanked_stacked,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											#xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
											xticks=(Array(1:num_years),labels),
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
		minx = 1 / num_games_total
		maxx = 1.0
		incx = (maxx - minx) / 5
		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
		incy = 1
		titlestring = L"\mbox{Number of teams eliminated over time}"
		xlabelstring = L"\mbox{Percent of season elapsed}"
		ylabelstring = L"\mbox{Number of teams eliminated}"
		legendtitlestring = "" #L"\mbox{Season}"
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
			#for s = 1:num_years
			#	curr_num_games = length(num_teams_eliminated[s,:])
			#	curr_label = labels[s]
      #  if length(col_labels) >= s
      #    plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), num_teams_eliminated[s,:], label=curr_label, color=col_labels[s]);
      #  else
      #    plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), num_teams_eliminated[s,:], label=curr_label);
      #  end
			#end
      curr_num_games = num_games_total
      curr_num_elim = sum(num_teams_eliminated, dims=1) / num_years
      errs = zeros(Float64, 2, num_games_total)
      min_num_elim = minimum(num_teams_eliminated, dims=1)
      max_num_elim = maximum(num_teams_eliminated, dims=1)
      errs[1,:] = curr_num_elim - min_num_elim
      errs[2,:] = max_num_elim - curr_num_elim
      curr_num_elim = curr_num_elim'
      plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), curr_num_elim, label="NBA average", color="black")
      plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), min_num_elim', color="black", linestyle="dashed", label="NBA min/max")
      plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), max_num_elim', color="black", linestyle="dashed")
      #errorbar(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), curr_num_elim, errs, color="black", linestyle="dashed", ecolor="black", errorevery=122)
			curr_num_games = length(avg_eliminated)
			plot(Array(1/curr_num_games:1/curr_num_games:curr_num_games/curr_num_games), avg_eliminated, label="simulated", color="blue", linestyle="dashdot") #marker="none", markevery=123)
			legend(loc="upper left", title=legendtitlestring)
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
		else
			fig = Plots.plot(show=false,
											title=titlestring,
											xlab=xlabelstring,
											ylab=ylabelstring,
											legendtitle=legendtitlestring,
											xticks=(Array(minx:incx:maxx),[@sprintf("\$%.0f\$", (100*i/num_games_total)) for i in minx:incx:maxx]),
											yticks=(Array(miny:incy:maxy),["\$$i\$" for i in miny:incy:maxy]),
											legend=:topleft,
											grid=false);
			for s = 1:size(num_teams_eliminated)[1]
				curr_num_games = length(num_teams_eliminated[s,:])
				curr_label = labels[s]
        if length(col_labels) >= s
          plot!(1:curr_num_games, num_teams_eliminated[s,:], label=curr_label, linecolor=col_labels[s]);
        else
          plot!(1:curr_num_games, num_teams_eliminated[s,:], label=curr_label);
        end
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

    writedlm(string(results_dir, "/", "nba_num_eliminated", csvext), num_teams_eliminated, ",")

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
	end # check if we should do_plotting

	return
end # main_parse

"""
    rankings_are_noisy

simulate a season to show that a finite set of samples is insufficient to get a good ranking estimate

NOTE: these results will match theory from closed_form_kendtau, but NOT from simulate code,
because simulate code uses a *tie-breaking* rule that when two teams have the same win percentage,
the head-to-head record is used to tie-break
"""
function rankings_are_noisy(;do_simulation=true, num_replications=1000, do_plotting=true, mode=MODE, results_dir="../results",
    num_teams = 30,
    num_playoff_teams = 16,
    num_rounds_set = [1 2 3 4 5 10 100 1000],
    # prob involves subdividing [0.5,1]
    prob = 0.5:0.5/50:1)
  Random.seed!(628) # for reproducibility

  ### DEBUG
  #num_rounds_set = [3]
  #prob = [.75]
  ##prob = [.5, .71375, .75]
  #num_teams = 3
  #num_playoff_teams = 0
  #set_mode(mode)
  ### DEBUG

	set_mode(mode)

	## Set constants
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2);
  num_prob = length(prob)

	## For output
	avg_kend = zeros(Float64, num_prob, length(num_rounds_set));

	if do_simulation
		for step_ind = 1:num_prob
			print("Step ", step_ind, " of ", num_prob, "\n")
			gamma = prob[step_ind];
			for num_rounds_ind in 1:length(num_rounds_set)
				num_rounds = num_rounds_set[num_rounds_ind];
				for rep = 1:num_replications
					stats = zeros(Int, num_teams, 2)
					for i = 1:num_teams
            stats[i,1] = i
						stats[i,2] = 0
					end

					for round_ind = 1:num_rounds
						for i = 1:num_teams
							for j = i+1:num_teams
								team_i_wins = teamWillWinNoTanking(i, j, gamma, true_strength, mode)

								# Do updates
								for k in [i,j]
									team_k_wins = (k == i) ? team_i_wins : !team_i_wins
									stats[k,2] += team_k_wins
								end
							end
						end
					end # loop over rounds

					## Calculate Kendell tau distance for this round
          curr_kend = kendtau(stats, 2, true_strength, mode, num_playoff_teams+1)
					avg_kend[step_ind, num_rounds_ind] += curr_kend / num_replications
          #println("wins = ", stats[:,2], "\tkend = $curr_kend\tavg_kend = $avg_kend")
				end # loop over replications
			end # loop over num_rounds_set
		end # loop over steps
    println("avg_kend = $avg_kend")
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
		incy = 10 #Int(floor((maxy - miny) / 5))
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

			close(fig)
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
end # rankings_are_noisy

"""
    model_validation
"""
function model_validation(;do_simulation = true, num_replications = 100000, 
    data_dir = "../data", results_dir = "../results", do_plotting = true,
    num_rounds = 3, num_steps = 2, gamma = 0.71425, 
    math_elim_mode = 0, selected_steps = nothing)
  Random.seed!(628) # for reproducibility
  selected_steps = clean_selected_steps(selected_steps)

  ## Simulation parameters
  #mode_list = [BT_ESTIMATED BT_DISTR]; mode_list_name = ["BT.est", "BT.beta"]
  mode_list = [BT_ESTIMATED]; mode_list_name= ["BT"]
  #mode_list = []; mode_list_name = []
  @assert(length(mode_list) == length(mode_list_name))
  #gamma_list = [0.50 0.55 0.60 0.65 0.70 0.7125 0.725 0.7375 0.75 0.80 0.85 0.90 0.95 1.00]
  #gamma_list = [0.7, 0.71, 0.715, 0.72, 0.75]
  #gamma_list = [0.7, 0.75]
  gamma_list = [0.7 0.71 0.711 0.7115 0.712 0.7125 0.713 0.7135 0.71375 0.714 0.71425 0.7145 0.715 0.72 0.725 0.75]
  num_modes = length(mode_list) + length(gamma_list)

  ## Stats we keep
  avg_stat    = 1
  stddev_stat = 2
  min_stat    = 3
  max_stat    = 4
  num_stats   = 4
  prefix      = ["avg_", "stddev_", "min_", "max_"]

  ## For output
  loss_list = zeros(Float64, num_modes, num_steps+1)
  win_pct_list = zeros(Float64, num_modes, num_steps+1, num_teams, num_stats)

  ## Data for comparison
  win_pct_nba = readdlm(string(data_dir, "/winpct.csv"), ',') # [year, team]
  num_header_rows = 1
  num_years = size(win_pct_nba, 2)

  for mode_ind = 1:num_modes
    curr_mode = mode_ind <= length(mode_list) ? mode_list[mode_ind] : STRICT
    curr_gamma = mode_ind <= length(mode_list) ? gamma : gamma_list[mode_ind - length(mode_list)]
    if curr_mode == STRICT
      println("Mode should be set to mode $mode_ind: ", curr_mode, " and gamma to $curr_gamma")
    else
      println("Mode should be set to mode $mode_ind: ", curr_mode)
    end
    set_mode(curr_mode, selected_steps)
    println("true_strength = ", true_strength)

    ## Retrieve win_pct matrix [step_ind, team_ind, stat]
    ## Save data
    if do_simulation
      win_pct = simulate(num_teams, num_playoff_teams, num_rounds, num_replications, num_steps, curr_gamma, breakpoint_list, nba_odds_list, nba_num_lottery, true_strength, curr_mode, math_elim_mode, selected_steps, true)
      println("win_pct = ", win_pct[:,:,avg_stat])
      win_pct_list[mode_ind, :, :, :] = win_pct
      for stat in [avg_stat]
        curr_name = "win_pct"
        if curr_mode == STRICT
          curr_name = string(curr_name, "_", curr_gamma)
        end
        writedlm(string(results_dir, "/", prefix[stat], curr_name, csvext), win_pct[:,:,stat], ',')
      end
    else
      win_pct = zeros(Float64, num_steps+1, num_teams, num_stats)
      for stat in [avg_stat]
        curr_name = "win_pct"
        if curr_mode == STRICT
          curr_name = string(curr_name, "_", curr_gamma)
        end
        win_pct[:,:,stat] = readdlm(string(results_dir, "/", prefix[stat], curr_name, csvext), ',')
        win_pct_list[mode_ind, :, :, stat] = win_pct[:,:,stat]
      end
    end # do_simulation

    ## Calculate mean squared error
    for step_ind = 1:num_steps+1
      loss = 0
      num_pts = num_teams * num_years
      for pos = 1:num_teams
        curr_calc = win_pct[step_ind, pos, avg_stat]
        for year = 1:num_years
          curr_real = win_pct_nba[num_header_rows + pos,year]
          loss += (curr_real-curr_calc)^2 / num_pts # mean squared error
        end # loop over years
      end # loop over teams
      loss_list[mode_ind, step_ind] = loss
      println("Step ", step_ind - 1, ": Loss from mode ", curr_mode, ": ", loss_list[mode_ind, step_ind])
    end # loop over steps
  end # iterate over modes in mode_list

  ## Save loss stats
  curr_name = "model_validation"
  writedlm(string(results_dir, "/", curr_name, csvext), loss_list, ',')
          
  ## Get avg nba data
  win_pct_nba_avg = sum(win_pct_nba[num_header_rows+1:num_header_rows+num_teams,:], dims=2)[:,1] / num_years
  errs = zeros(Float64, 2, num_teams)
  errs[2,:] = maximum(win_pct_nba[num_header_rows+1:num_header_rows+num_teams,:], dims=2) - win_pct_nba_avg
  errs[1,:] = win_pct_nba_avg - minimum(win_pct_nba[num_header_rows+1:num_header_rows+num_teams,:], dims=2)

  ## Plot simulated vs real average win pct, with error bars
  if do_plotting
    print("Plotting win_pct\n")
    minx = 1
    incx = 5
    maxx = num_teams
    miny = 0
    incy = 10
    maxy = 100
    titlestring = L"\mbox{Win percentage by rank}"
    xlabelstring = L"\mbox{Rank of team}"
    ylabelstring = L"\mbox{Winning percentage at end of season}"
    legendtitlestring = L"\mbox{Method}"
    fname_stub = "win_pct"
    if use_pyplot
      for tank_ind in 1:num_steps+1 #[1,num_steps+1]
        tank_name = string("_",tank_ind-1,"tank")
        fname = string(results_dir,"/",ext_folder,"/",fname_stub,tank_name,ext)
        fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,tank_name,lowext)
        fig = figure(frameon=false)
        title(titlestring)
        xlabel(xlabelstring)
        ylabel(ylabelstring)
        #xticks(Array(minx:incx:maxx))
        tmp = Array(minx-1:incx:maxx)
        tmp[1] += 1
        xticks(tmp)
        yticks(Array(miny:incy:maxy))
        gammas_to_plot = [0.71425]
        gammas_to_plot_ind = zeros(Int, length(gammas_to_plot))
        for r = 1:length(gammas_to_plot)
          tmp = findfirst(isequal(gammas_to_plot[r]), gamma_list)
          if !isa(tmp, Nothing)
            gammas_to_plot_ind[r] = tmp[2]
          end
        end
        num_modes_to_plot = length(mode_list) + length(gammas_to_plot)
        for r = 1:num_modes_to_plot
          if r <= length(mode_list)
            curr_style="dashdot"
            curr_marker=""
            curr_size=0
          elseif r - length(mode_list) <= length(style)
            tmp_ind = r - length(mode_list)
            curr_style=style[tmp_ind]
            curr_marker=shape[tmp_ind]
            curr_size=shapesize[tmp_ind]
          else
            curr_style="solid"
            curr_marker=""
            curr_size=0
          end
          curr_ind = (r <= length(mode_list)) ? r : gammas_to_plot_ind[r-length(mode_list)]
          if curr_ind <= 0
            continue
          end
          curr_label = (r <= length(mode_list)) ? mode_list_name[r] : latexstring("\\gamma=",gamma_list[curr_ind])
          plot(1:num_teams, win_pct_list[curr_ind,tank_ind,:,avg_stat], label=curr_label, linestyle=curr_style, marker=curr_marker, markersize=curr_size)
        end
        curr_label = "NBA average"
        plot(1:num_teams, win_pct_nba_avg, label=curr_label, color="black", marker="", markersize=5)
        errorbar(1:num_teams, win_pct_nba_avg, errs, color="black")

        legend(loc="upper right", title=legendtitlestring) 
        #legend(bbox_to_anchor=[0.5,.95, 0.5, 0.45], title=legendtitlestring, ncol=4, fontsize="xx-small", markerscale=0.4)
        #legend(loc="upper right", title=legendtitlestring, ncol=4, fontsize="xx-small", markerscale=0.4)
        PyPlot.savefig(fname)
        PyPlot.savefig(fname_low)
        close(fig)
      end # loop over tank indices we want to plot
    end # check if pyplot is used

    ## Plot loss
    println("Plotting loss")
    miny = Int(floor(minimum(loss_list)))
    incy = 1
    maxy = Int(ceil(maximum(loss_list)))
		titlestring = L"\mbox{Model error for Bradley-Terry MLE and various choices of $\gamma$}"
    xlabelstring = L"\mbox{Model (Bradley-Terry or value of $\gamma$)}"
		ylabelstring = L"\mbox{Mean squared error from NBA data}"
		#legendtitlestring = L"\mbox{Number of selfish teams}"
    legendtitlestring = ""
		fname_stub = "model_loss"
		fname = string(results_dir,"/",ext_folder,"/",fname_stub,ext)
		fname_low = string(results_dir,"/",lowext_folder,"/",fname_stub,lowext)

    if use_pyplot
			fig = figure(frameon=false)
			title(titlestring)
      xlabel(xlabelstring)
			ylabel(ylabelstring)
			#xticks(Array(minx:incx:maxx))
      gamma_list_name = [string(gamma_list[i]) for i = 1:length(gamma_list)]
      xticks(1:num_modes, vcat(mode_list_name, gamma_list_name),rotation=-30)
			yticks(Array(miny:incy:maxy))
      for r = 1:num_steps+1
        curr_num = Int(num_teams * (r-1) / num_steps)
        curr_label = string("$curr_num selfish teams")
        #plot(1:num_modes, loss_list[:,r], label=curr_label, color=col[r], linestyle=style[r], marker="", markersize=5)
        plot(1:num_modes, loss_list[:,r], label=curr_label, color=col[r], linestyle="none", marker=shape[r], markersize=shapesize[r])
      end
      legend(loc="upper center", title=legendtitlestring) 
			PyPlot.savefig(fname)
			PyPlot.savefig(fname_low)
			close(fig)
    end # check if pyplot is used
  end # check if do_plotting
    
  println("## Summary of model validation experiments ##")
  for mode_ind = 1:num_modes
    curr_mode = mode_ind <= length(mode_list) ? mode_list[mode_ind] : STRICT
    curr_gamma = mode_ind <= length(mode_list) ? gamma : gamma_list[mode_ind - length(mode_list)]
    for step_ind = 1:num_steps+1 
      if curr_mode == STRICT
        println("Step ", step_ind - 1, ": Loss from mode ", curr_mode, " (gamma: ", curr_gamma, "): ", loss_list[mode_ind, step_ind])
      else
        println("Step ", step_ind - 1, ": Loss from mode ", curr_mode, ": ", loss_list[mode_ind, step_ind])
      end
    end
  end
end # model_validation

"""
    closed_form_kendtau
"""
function closed_form_kendtau(;
    num_teams = 30,
    num_playoff_teams = 16,
    gamma = 0.71425,
    num_rounds = 3,
    mode = MODE)
  Random.seed!(628) # for reproducibility
  set_mode(mode)
  eps_diff = 1/num_teams
  
  num_pairs = Int(num_teams * (num_teams-1) / 2)
  num_perm = (num_rounds + 1)^num_pairs
  @assert(num_perm < 1e10 && num_perm > 0) # make sure we do not have too many possibilities
  # Win totals can increment by (0,num_rounds),(1,num_rounds-1),...,(num_rounds,0) for each pair of teams that plays
  pair_wins = zeros(Int, num_pairs)
  pair_ind = 1
  avg_kend = 0
  for perm = 1:num_perm
    prob = 1
    stats = Matrix{Any}(undef, num_teams, 2)
    for i = 1:num_teams
      stats[i,1] = i
      stats[i,2] = 0.0
    end

    tmp_ind = 1
    for i = 1:num_teams
      for j = i+1:num_teams
        wins_i = pair_wins[tmp_ind]
        wins_j = num_rounds - wins_i
        stats[i,2] += wins_i
        stats[j,2] += wins_j
        prob *= gamma^wins_i * (1-gamma)^wins_j * binomial(num_rounds, wins_i)
        tmp_ind += 1
      end
    end
    
    # Now we calculated Kendall tau distance when the end-of-season win totals are as given
    # The only thing we have to be careful about is ties among teams
    # We should consider all possible permutations of the tied teams 
    # To do this, we put together all equivalence classes, where these are based on number of wins each team has
    # (range is from 0 to num_rounds * (num_teams-1))
    max_num_wins = 1 + num_rounds * (num_teams-1)
    equiv_classes = Array{Array{Int64,1}}(undef, max_num_wins)
    for num_wins = 1:max_num_wins
      equiv_classes[num_wins] = []
    end
    num_equiv_perm = 1
    for i = 1:num_teams
      curr_wins = Int(stats[i,2])
      equiv_classes[curr_wins+1] = [equiv_classes[curr_wins+1]; i]
      num_equiv_perm *= length(equiv_classes[curr_wins+1])
    end
    equiv_class_perm = Array{Any}(undef, max_num_wins)
    for num_wins = 1:max_num_wins
      curr_perms = collect(permutations(equiv_classes[num_wins]))
      equiv_class_perm[num_wins] = curr_perms
    end
    equiv_perm_ind = ones(Int64, max_num_wins)
    for tmp_perm = 1:num_equiv_perm
      for tmp_ind = 1:max_num_wins
        num_teams_in_class = length(equiv_classes[tmp_ind])
        curr_perm = equiv_perm_ind[tmp_ind]
        for tmp_team = 1:num_teams_in_class
          curr_team = equiv_class_perm[tmp_ind][curr_perm][tmp_team]
          stats[curr_team, 2] = floor(stats[curr_team,2]) + eps_diff * (tmp_team - 1) / num_teams_in_class
        end
      end
      curr_kend = kendtau(stats, 2, true_strength, mode, num_playoff_teams+1)
      avg_kend += prob * curr_kend / num_equiv_perm

      ### DEBUG
      #println("wins = ", stats[:,2])
      #println("prob = $prob") 
      #println("kend = ", curr_kend)

      # Go to the next permutation
      for tmp_ind = 1:max_num_wins
        if equiv_perm_ind[tmp_ind] < length(equiv_class_perm[tmp_ind])
          equiv_perm_ind[tmp_ind] += 1
          break
        else
          equiv_perm_ind[tmp_ind] = 1
        end
      end
    end
    #println("perm $perm, avg_kend = ", avg_kend)

    # Increment pair_wins
    continue_inc = true
    while continue_inc
      if pair_wins[pair_ind] < num_rounds
        pair_wins[pair_ind] += 1
        pair_ind = 1
        continue_inc = false
      else
        pair_wins[pair_ind] = 0
        pair_ind += 1
        if pair_ind > num_pairs
          continue_inc = false
        end
      end
    end # while continue_inc
    if pair_ind > num_pairs
      break
    end
  end # loop over perm
  println("avg_kend = $avg_kend")
end # closed_form_kendtau

"""
    count_num_win_partitions

Count the number of win totals `n` teams can have when each team can have at most `M` wins
and the sum of all wins must equal `total`
"""
function count_num_win_partitions(total, n, M)
  if (M * n < total) || (n <= 0) || (M < 0) || (total < 0)
    return 0
  elseif (M * n == total) || (n == 1)
    return 1
  else
    tmp_sum = 0
    for i = 0:M
      tmp_sum += count_num_win_partitions(total - i, n-1, M)
    end
    return tmp_sum
  end
end # count_num_win_partitions

"""
    tanking_unit_tests
"""
function tanking_unit_tests()
	test_ranking = 1:num_teams
	test_strength = 1:num_teams
	@assert ( kendtau_sorted(test_ranking, test_strength, 1) == Int(num_teams * (num_teams-1) / 2) )
	test_strength = num_teams:-1:1
	@assert ( kendtau_sorted(test_ranking, test_strength, 1) == 0 )
end # tanking_unit_tests

end # module Tanking
