# Needs Plots, DelimitedFiles, GR
using Plots
using StatPlots
using DelimitedFiles
using LaTeXStrings
using Printf
include("simulate.jl")
include("parse.jl")

## When the (draft) ranking will be set as fraction of number games
set_ranking = [4//8; 5//8; 6//8; 7//8; 1]
num_rankings = length(set_ranking)
shape = [:vline, :utriangle, :rect, :x, :circle]
col = ["red", "orange", "green", "blue", "black"]
#color_for_cutoff_point = ["c" "b" "m" "r" "k"]

## For plotting
ext = ".svg"

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
		## Do plotting
		#pgfplots() # Pkg.add("PGFPlots")
		#pyplot() # Pkg.add("PyPlot") Pkg.add("PyCall") Pkg.add("LaTeXStrings")
		gr(dpi=200); # Pkg.add("GR")
		#gr(); # Pkg.add("GR")
		upscale = 1 # upscaling in resolution
		fntsm = Plots.font("sans-serif", 8.0 * upscale)
		fntlg = Plots.font("sans-serif", 12.0 * upscale)
		#default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
		default(size=(600*upscale,400*upscale)) # plot canvas size

		## Plot avg_kend (Kendell tau distance)
		print("Plotting avg_kend: average swap distance\n")
		miny = Int(ceil(findmin(avg_kend)[1]))
		maxy = Int(floor(findmax(avg_kend)[1]))

		fig = Plots.plot(show=false,
										#title=L"\mbox{Fidelity of ranking by breapoint and tanking probability}",
										title=L"\mbox{Accuracy of ranking}",
										#xlab=L"\mbox{Percentage of tanking teams}", 
										xlab=L"\mbox{Probability of tanking}", 
										#ylab=latexstring("Test1","\\mbox{}\\\\","Test2"),
										#ylab=latexstring("\\mbox{Kendell tau distance from}", "\\\\", "\\mbox{true ranking of non-playoff teams}"),
										ylab=L"\mbox{Distance from true ranking of non-playoff teams}",
										#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
										xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
										yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
										legend=:bottomright,
										#legend=:best,
										legendfont=6,
										legendtitle=L"\mbox{Draft ranking breakpoint}",
										titlefont=12,
										tickfont=8,
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
		savefig(fig, string(results_dir,"/avg_kend",ext));

		## Plot avg_games_tanked (# games tanked by draft ranking breakpoint)
		print("Plotting avg_games_tanked: average number of games tanked\n")
		miny = Int(ceil(findmin(avg_games_tanked)[1]))
		maxy = Int(floor(findmax(avg_games_tanked)[1]))
		inc = (maxy - miny) / 5
		fig = Plots.plot(show=false,
										title=L"\mbox{Total games tanked}",
										#xlab=L"\mbox{Percentage of tanking teams}", 
										xlab=L"\mbox{Probability of tanking}", 
										ylab=L"\mbox{Number of tanked games}",
										#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
										xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
										yticks=(Array(miny:inc:maxy),[@sprintf("\$%.0f\$", i) for i in miny:inc:maxy]),
										legend=:topleft,
										#legend=:best,
										legendfont=6,
										legendtitle=L"\mbox{Draft ranking breakpoint}",
										titlefont=12,
										tickfont=8,
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
		savefig(fig, string(results_dir,"/avg_games_tanked",ext));

		## Plot avg_already_tank
		print("Plotting avg_already_tank: average number of tanking teams\n")
		miny = Int(ceil(findmin(avg_already_tank)[1]))
		maxy = Int(floor(findmax(avg_already_tank)[1]))
		fig = Plots.plot(show=false,
										title=L"\mbox{Number of tanking teams}",
										#xlab=L"\mbox{Percentage of tanking teams}", 
										xlab=L"\mbox{Probability of tanking}", 
										ylab=L"\mbox{Average number of tanking teams}",
										#xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
										xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
										yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
										legend=:topleft,
										legendfont=6,
										legendtitle=L"\mbox{Draft ranking breakpoint}",
										titlefont=12,
										tickfont=8,
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
		savefig(fig, string(results_dir,"/avg_already_tank",ext));

		## Plot avg_eliminated
		print("Plotting avg_eliminated: average number of eliminated teams by every game of the season\n")
		miny = Int(ceil(findmin(avg_eliminated)[1]))
		maxy = Int(floor(findmax(avg_eliminated)[1]))
		inc = num_games / 5
		fig = Plots.plot(show=false,
										title=L"\mbox{Number teams eliminated over time}",
										xlab=L"\mbox{Percent through season}", 
										ylab=L"\mbox{Number teams eliminated}",
										xticks=(Array(0:inc:num_games),[@sprintf("\$%.0f\$", (100*i/num_games)) for i in 0:inc:num_games]),
										yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
										#legend=:topleft,
										#legendfont=6,
										#legendtitle=L"\mbox{Tanking probability}",
										legend=:none,
										titlefont=12,
										tickfont=8,
										grid=false);
			#step_ind = Int(floor((num_steps+1) / 2))
			#tank_perc = 100 * (step_ind - 1) / num_steps
			#plot!(1:num_games, avg_eliminated[step_ind,:], label="\$$tank_perc\$");
		plot!(1:num_games, sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1), linetype=:bar);
		#plot!(1:num_games, sum(avg_eliminated, dims=1)[1,:] / (num_steps + 1));
		savefig(fig, string(results_dir,"/avg_eliminated",ext));
	end # if do_plotting
	return
end # main_simulate

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
		ind = [3,4,5]
		num_years=5

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

		num_games_tanked_stacked = zeros(Int, num_years, length(set_ranking))
		num_games_tanked_stacked[:,1] = num_games_tanked[:,1]
		num_games_tanked_stacked[:,2] = num_games_tanked[:,2] - num_games_tanked[:,1]
		num_games_tanked_stacked[:,3] = num_games_tanked[:,3] - num_games_tanked[:,2]
		num_games_tanked_stacked[:,4] = num_games_tanked[:,4] - num_games_tanked[:,3]
		num_games_tanked_stacked[:,5] = num_games_tanked[:,5] - num_games_tanked[:,4]

		miny = Int(ceil(findmin(num_games_tanked)[1]))
		maxy = Int(floor(findmax(num_games_tanked)[1]))
		inc = (maxy - miny) / 5
		ctg = repeat(vcat([latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}") for r in ind if r < length(set_ranking)], L"\mbox{end of season}"), inner=num_years)
		fig = groupedbar(num_games_tanked_stacked[:,ind], 
										xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
										yticks=(Array(0:50:maxy),[@sprintf("\$%d\$", i) for i in 0:50:maxy]),
										group=ctg,
										lw=0,
										#bar_position=:dodge,
										bar_position=:stack,
										title=L"\mbox{Number of games that could be tanked}",
										xlab=L"\mbox{Season}", 
										ylab=L"\mbox{Number of possibly tanked games}",
										legend=:best,
										legendfont=6,
										legendtitle=L"\mbox{Draft ranking breakpoint}",
										titlefont=12,
										tickfont=8,
										grid=false,
										show=false)
		savefig(fig, string(results_dir,"/nba_num_games_tanked",ext));

		# Print number of teams eliminated
		print("Plotting num_teams eliminated: number of eliminated teams by the breakpoint mark\n")
		#num_teams_eliminated = zeros(Int, num_years, length(set_ranking))
		#num_teams_eliminated[1,:] = num_teams_eliminated_1314
		#num_teams_eliminated[2,:] = num_teams_eliminated_1415
		#num_teams_eliminated[3,:] = num_teams_eliminated_1516
		#num_teams_eliminated[4,:] = num_teams_eliminated_1617
		#num_teams_eliminated[5,:] = num_teams_eliminated_1718
		num_teams_eliminated = hcat(num_teams_eliminated_1314, num_teams_eliminated_1415, num_teams_eliminated_1516, num_teams_eliminated_1617, num_teams_eliminated_1718)
		num_teams_eliminated = num_teams_eliminated'
		#print(num_teams_eliminated,"\n")

		num_teams_eliminated_stacked = zeros(Int, num_years, length(set_ranking))
		num_teams_eliminated_stacked[:,1] = num_teams_eliminated[:,1]
		num_teams_eliminated_stacked[:,2] = num_teams_eliminated[:,2] - num_teams_eliminated[:,1]
		num_teams_eliminated_stacked[:,3] = num_teams_eliminated[:,3] - num_teams_eliminated[:,2]
		num_teams_eliminated_stacked[:,4] = num_teams_eliminated[:,4] - num_teams_eliminated[:,3]
		num_teams_eliminated_stacked[:,5] = num_teams_eliminated[:,5] - num_teams_eliminated[:,4]

		miny = Int(ceil(findmin(num_teams_eliminated)[1]))
		maxy = Int(floor(findmax(num_teams_eliminated)[1]))
		inc = (maxy - miny) / 5
		ctg = repeat(vcat([latexstring(numerator(set_ranking[r]),"/",denominator(set_ranking[r]), "\\mbox{ through season}") for r in ind if r < length(set_ranking)], L"\mbox{end of season}"), inner=num_years)
		fig = groupedbar(num_teams_eliminated[:,ind], 
										xticks=(Array(1:num_years),["\$13-14\$","\$14-15\$","\$15-16\$","\$16-17\$","\$17-18\$"]),
										#yticks=(Array(miny:inc:maxy),[@sprintf("\$%d\$", i) for i in miny:inc:maxy]),
										group=ctg,
										lw=0,
										#bar_position=:dodge,
										bar_position=:stack,
										title=L"\mbox{Number of teams eliminated by breakpoint mark}",
										xlab=L"\mbox{Season}", 
										ylab=L"\mbox{Number of teams eliminated}",
										legend=:best,
										legendfont=6,
										legendtitle=L"\mbox{Draft ranking breakpoint}",
										titlefont=12,
										tickfont=8,
										grid=false,
										show=false)
		savefig(fig, string(results_dir,"/nba_num_teams_eliminated",ext));
	end

	return
end # main_parse
