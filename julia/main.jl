## Required dependencies
using DelimitedFiles
using LaTeXStrings
using Printf
include("simulate.jl")
include("parse.jl")

TITLE_FONTSIZE=10
AXIS_TITLE_FONTSIZE=10
TICK_LABEL_FONTSIZE=8
LEGEND_FONTSIZE=8
LEGEND_TITLE_FONTSIZE=8
DPI=200

## For plotting
use_pyplot = true
ext_folder = "pdf"
lowext_folder = "png"
ext = string(".",ext_folder)
lowext = string("_low.",lowext_folder)
upscale = 1 # upscaling in resolution
if !use_pyplot
	#ext = ".svg"
	using Plots
	using StatPlots
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
	using PyCall
	pygui(:qt5) # others do not work on mac
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


## When the (draft) ranking will be set as fraction of number games
#set_ranking = [4//8; 5//8; 6//8; 7//8; 1]
set_ranking = [1//2; 2//3; 3//4; 5//6; 7//8; 1]
num_rankings = length(set_ranking)
shape = [:vline, :utriangle, :rect, :x, :triangle, :circle]
col = ["red", "orange", "green", "blue", "violet", "black"]
#color_for_cutoff_point = ["c" "b" "m" "r" "k"]

function main_simulate(do_simulation = true, num_repeats = 100000, do_plotting=true, results_dir = "../results")
	## Variables that need to be set
	num_teams = 30 # number of teams
	num_rounds = 3 # a round consists of each team playing each other team
	num_steps = 20 # discretization of [0,1] for tanking probability
	gamma = 0.75 # probability a better-ranked team wins over a worse-ranked team
	## end variables that need to be set

	## Set constants
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2)
	num_games = num_rounds * num_games_per_round

	## For output
	avg_kend = 0 # [step,cutoff], holds KT distance for each tanking probability and cutoff for draft ranking
	avg_games_tanked = 0 # [step,cutoff], number games tanked by cutoff
	avg_already_tank = 0 # [step,cutoff], number teams already tanking by the cutoff
	avg_eliminated = 0 # [step,game], number teams eliminated by each game

	## Do simulation or retrieve data
	if do_simulation
		## Do simulation
		avg_kend, avg_games_tanked, avg_already_tank, avg_eliminated = simulate(num_teams, num_rounds, num_repeats, num_steps, gamma, set_ranking)
		writedlm(string(results_dir, "/avg_kend.csv"), avg_kend, ',')
		writedlm(string(results_dir, "/avg_games_tanked.csv"), avg_games_tanked, ',')
		writedlm(string(results_dir, "/avg_already_tank.csv"), avg_already_tank, ',')
		writedlm(string(results_dir, "/avg_eliminated.csv"), avg_eliminated, ',')
	else
		avg_kend = readdlm(string(results_dir, "/avg_kend.csv"), ',')
		avg_games_tanked = readdlm(string(results_dir, "/avg_games_tanked.csv"), ',')
		avg_already_tank = readdlm(string(results_dir, "/avg_already_tank.csv"), ',')
		avg_eliminated = readdlm(string(results_dir, "/avg_eliminated.csv"), ',')
		num_steps = size(avg_kend)[1] - 1
	end

	if (do_plotting)
		## Plot avg_kend (Kendell tau distance)
		print("Plotting avg_kend: average swap distance\n")
		minx = 0
		incx = 0.1
		maxx = 1
		miny = Int(floor(findmin(avg_kend)[1]))
		incy = 1
		maxy = Int(ceil(findmax(avg_kend)[1]))
		#titlestring=L"\mbox{Fidelity of ranking by breapoint and tanking probability}"
		titlestring = L"\mbox{Accuracy of ranking}"
		xlabelstring = L"\mbox{Probability of tanking once eliminated}"
		ylabelstring = L"\mbox{Distance from true ranking of non-playoff teams}"
		legendtitlestring = L"\mbox{Draft ranking breakpoint}"
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
				if set_ranking[r] == 1	
					curr_label = "end of season" #L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot(0:(1/num_steps):1, avg_kend[:,r], label=curr_label, color=col[r])
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
				#cutoff_game = set_ranking[r]
				curr_label = ""
				if set_ranking[r] == 1	
					#curr_label = latexstring("$tmp", "\\mbox{ through season}")
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot!(0:(1/num_steps):1, avg_kend[:,r], label=curr_label, linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot avg_games_tanked (# games tanked by draft ranking breakpoint)
		print("Plotting avg_games_tanked: average number of games tanked\n")
		miny = Int(floor(findmin(avg_games_tanked)[1]))
		incy = 50 #Int(ceil((maxy - miny) / (5 * 10)) * 10)
		maxy = Int(ceil(findmax(avg_games_tanked)[1] / incy) * incy)  #Int(floor(findmax(avg_games_tanked)[1]))
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
				if set_ranking[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot(0:(1/num_steps):1, avg_games_tanked[:,r], label=curr_label, color=col[r])
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
				if set_ranking[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot!(0:(1/num_steps):1, avg_games_tanked[:,r], label=curr_label, linecolor=col[r]);
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
		miny = Int(floor(findmin(avg_already_tank)[1]))
		incy = 1
		maxy = Int(ceil(findmax(avg_already_tank)[1]))
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
				if set_ranking[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot(0:(1/num_steps):1, avg_already_tank[:,r], label=curr_label, color=col[r])
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
				if set_ranking[r] == 1	
					curr_label = L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
				end
				plot!(0:(1/num_steps):1, avg_already_tank[:,r], label=curr_label, linecolor=col[r]);
								#markershape=shape[r], markersize=2, markercolor=col[r], markerstrokecolor=col[r]);
			end
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end

		## Plot avg_eliminated
		print("Plotting avg_eliminated: average number of effectively eliminated teams by every game of the season\n")
		minx = 1
		maxx = num_games
		incx = (maxx - minx) / 5
		miny = Int(floor(findmin(avg_eliminated)[1]))
		incy = 1
		maxy = Int(ceil(findmax(avg_eliminated)[1]))
		titlestring = L"\mbox{Number teams effectively eliminated over time}"
		xlabelstring = L"\mbox{Percent through season}"
		ylabelstring = L"\mbox{Number teams effectively eliminated}"
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
			y = sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1)
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
			plot!(1:num_games, sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1), linetype=:bar);
			Plots.savefig(fname)
			Plots.savefig(fname_low)
		end
	end # if do_plotting
	return
end; # main_simulate

function main_parse(do_plotting=true, data_dir="../data", results_dir="../results")
	num_teams_eliminated_1314, num_games_tanked_1314, stats1314, critical_game1314 = parseNBASeason("games1314.csv", set_ranking, data_dir)
	num_teams_eliminated_1415, num_games_tanked_1415, stats1415, critical_game1415 = parseNBASeason("games1415.csv", set_ranking, data_dir)
	num_teams_eliminated_1516, num_games_tanked_1516, stats1516, critical_game1516 = parseNBASeason("games1516.csv", set_ranking, data_dir)
	num_teams_eliminated_1617, num_games_tanked_1617, stats1617, critical_game1617 = parseNBASeason("games1617.csv", set_ranking, data_dir)
	num_teams_eliminated_1718, num_games_tanked_1718, stats1718, critical_game1718 = parseNBASeason("games1718.csv", set_ranking, data_dir)

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

	if (do_plotting)
		ind = [3,5,6] # needs to be ascending
		@assert ( length(set_ranking) in ind )
		num_years=5
		labels = [L"2013-2014", L"2014-2015", L"2015-2016", L"2016-2017", L"2017-2018"]
		col_labels = ["red", "orange", "green", "blue", "black"]

		## Plot # games tanked
		print("Plotting num_games_tanked: number of games (possibly) tanked by the breakpoint mark\n")
		#num_games_tanked = zeros(Int, num_years, length(set_ranking))
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
		legendtitlestring = L"\mbox{Draft ranking breakpoint}"
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
				if set_ranking[r] == 1
					curr_label=L"\mbox{end of season}"
				else
					curr_label = latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}")
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
			ctg = repeat(vcat([latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}") for r in ind if r < length(set_ranking)], L"\mbox{end of season}"), inner=num_years)
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
		print("Plotting num_teams_eliminated: number of effectively eliminated teams by every game of the season\n")
		num_games = 15 * 82;
		num_teams_eliminated = hcat(num_teams_eliminated_1314, num_teams_eliminated_1415, num_teams_eliminated_1516, num_teams_eliminated_1617, num_teams_eliminated_1718)
		num_teams_eliminated = num_teams_eliminated'

		minx = 1
		maxx = num_games
		incx = (maxx - minx) / 5
		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
		incy = 1
		titlestring = L"\mbox{Number teams effectively eliminated over time}"
		xlabelstring = L"\mbox{Percent through season}"
		ylabelstring = L"\mbox{Number teams effectively eliminated}"
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
			yticks(Array(miny:incy:maxy))
			for s = 1:size(num_teams_eliminated)[1]
				curr_num_games = length(num_teams_eliminated[s,:])
				curr_label = labels[s]
				plot(1:curr_num_games, num_teams_eliminated[s,:], label=curr_label, color=col_labels[s]);
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
#		#num_teams_eliminated = zeros(Int, num_years, length(set_ranking))
#		#num_teams_eliminated[1,:] = num_teams_eliminated_1314
#		#num_teams_eliminated[2,:] = num_teams_eliminated_1415
#		#num_teams_eliminated[3,:] = num_teams_eliminated_1516
#		#num_teams_eliminated[4,:] = num_teams_eliminated_1617
#		#num_teams_eliminated[5,:] = num_teams_eliminated_1718
#		num_teams_eliminated = hcat(num_teams_eliminated_1314, num_teams_eliminated_1415, num_teams_eliminated_1516, num_teams_eliminated_1617, num_teams_eliminated_1718)
#		num_teams_eliminated = num_teams_eliminated'
#		#print(num_teams_eliminated,"\n")
#
#		num_teams_eliminated_stacked = zeros(Int, num_years, length(set_ranking))
#		num_teams_eliminated_stacked[:,1] = num_teams_eliminated[:,1]
#		num_teams_eliminated_stacked[:,2] = num_teams_eliminated[:,2] - num_teams_eliminated[:,1]
#		num_teams_eliminated_stacked[:,3] = num_teams_eliminated[:,3] - num_teams_eliminated[:,2]
#		num_teams_eliminated_stacked[:,4] = num_teams_eliminated[:,4] - num_teams_eliminated[:,3]
#		num_teams_eliminated_stacked[:,5] = num_teams_eliminated[:,5] - num_teams_eliminated[:,4]
#
#		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
#		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
#		inc = floor((maxy - miny) / 5)
#		ctg = repeat(vcat([latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}") for r in ind if r < length(set_ranking)], L"\mbox{end of season}"), inner=num_years)
#		fig = groupedbar(num_teams_eliminated[:,ind], 
#										xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
#										yticks=(Array(miny:inc:maxy),[@sprintf("\$%d\$", i) for i in miny:inc:maxy]),
#										group=ctg,
#										lw=0,
#										bar_position=:dodge,
#										#bar_position=:stack,
#										title=L"\mbox{Number of teams effectively eliminated by breakpoint mark}",
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

function rankings_are_noisy(do_simulation=true, results_dir="../results")
	num_steps = 50; # number of subdivisions of [0.5,1]
	prob = 0.5:0.5/num_steps:1;
	num_teams = 30;
	num_repeats = 1000;
	num_rounds_set = [1 2 3 4 5 10 100 1000];
	num_games_per_round = Int(num_teams * (num_teams - 1) / 2);

	avg_kend = zeros(Float64, num_steps+1, length(num_rounds_set));

	if do_simulation
		for step_ind = 1:num_steps+1
			print("Step ", step_ind, " of ", num_steps+1, "\n")
			gamma = prob[step_ind];
			for num_rounds_ind in 1:length(num_rounds_set)
				num_rounds = num_rounds_set[num_rounds_ind];
				for rep = 1:num_repeats
					stats = zeros(Int, num_teams, 2);
					for i = 1:num_teams
						stats[i,1] = i;
						stats[i,2] = 0;
					end

					for round_ind = 1:num_rounds
						for i = 1:num_teams
							for j = i+1:num_teams
								if rand() < gamma
									team_i_wins = true
								else
									team_i_wins = false
								end

								# Do updates
								for k in [i,j]
									team_k_wins = (k == i) ? team_i_wins : !team_i_wins
									stats[k,2] = stats[k,2] + team_k_wins
								end
							end
						end
					end # loop over rounds

					## Calculate Kendell tau distance for this round
					avg_kend[step_ind, num_rounds_ind] = avg_kend[step_ind, num_rounds_ind] + kendtau(stats,2) / num_repeats
				end # loop over repeats
			end # loop over num_rounds_set
		end # loop over steps
		writedlm(string(results_dir, "/noisy_ranking.csv"), avg_kend, ',')
	else
		avg_kend = readdlm(string(results_dir, "/noisy_ranking.csv"), ',')
	end

  ## Plot noisy ranking
	print("Plotting noisy_ranking: shows low fidelity of ranking unless teams play each other many times\n")
	minx = 0.5
	incx = 0.1
	maxx = 1
	miny = Int(ceil(findmin(avg_kend)[1]));
	incy = 50 #Int(floor((maxy - miny) / 5))
	maxy = Int(floor(findmax(avg_kend)[1]));
	#maxy = Int(ceil(maxy/incy)*incy)
	titlestring = L"\mbox{Accuracy of ranking}"
	xlabelstring = L"\mbox{Probability better team wins}"
	ylabelstring = L"\mbox{Distance from true ranking}"
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
end; # rankings_are_noisy
