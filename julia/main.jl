# Needs Plots, DelimitedFiles, GR
using Plots
using DelimitedFiles
using LaTeXStrings
using Printf
include("tanking.jl")

function main(do_simulation = true, num_repeats = 100000, results_dir = "../results")
	## Variables that need to be set
	num_teams = 30 # number of teams
	num_rounds = 3 # a round consists of each team playing each other team
	num_steps = 20 # discretization of [0,1] for tanking probability
	gamma = 0.75 # probability a better-ranked team wins over a worse-ranked team

	## When the (draft) ranking will be set as fraction of number games
	set_ranking = [4//8; 5//8; 6//8; 7//8; 1]
	num_rankings = length(set_ranking)
	col = ["red", "orange", "green", "blue", "black"]
	shape = [:vline, :utriangle, :rect, :x, :circle]
	#color_for_cutoff_point = ["c" "b" "m" "r" "k"]
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

	## Do plotting
	#pgfplots() # Pkg.add("PGFPlots")
	#pyplot() # Pkg.add("PyPlot") Pkg.add("PyCall") Pkg.add("LaTeXStrings")
	#gr(dpi=300); # Pkg.add("GR")
	gr(); # Pkg.add("GR")
	upscale = 1 # upscaling in resolution
	fntsm = Plots.font("sans-serif", 8.0 * upscale)
	fntlg = Plots.font("sans-serif", 12.0 * upscale)
	#default(titlefont=fntlg, guidefont=fntlg, tickfont=fntsm, legendfont=fntsm)
	#default(size=(800*upscale,600*upscale)) # plot canvas size

	## Plot avg_kend (Kendell tau distance)
	print("Plotting avg_kend: average swap distance\n")
	miny = Int(ceil(findmin(avg_kend)[1]))
	maxy = Int(floor(findmax(avg_kend)[1]))

	fig = Plots.plot(show=false,
			xlab=L"\mbox{Probability of tanking once eliminated}", 
			ylab=L"\mbox{Distance from true ranking of non-playoff teams}",
			title=L"\mbox{Kendall tau distance by cutoff and tanking probability}",
			xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
			yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
			legend=:bottomright,
			#legend=:best,
			legendfont=6,
			legendtitle=L"\mbox{Draft ranking cutoff}",
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
	savefig(fig, string(results_dir,"/avg_kend",".pdf"));
	#savefig(fig, "plot.pdf");
	#save(string(results_dir,"/plot",".pdf"), fig);

	## Plot avg_games_tanked (# games tanked by draft ranking cutoff)
	print("Plotting avg_games_tanked: average number of games tanked\n")
	miny = Int(ceil(findmin(avg_games_tanked)[1]))
	maxy = Int(floor(findmax(avg_games_tanked)[1]))
	inc = (maxy - miny) / 5
	fig = Plots.plot(show=false,
			xlab=L"\mbox{Probability of tanking once eliminated}", 
			ylab=L"\mbox{Number of tanked games}",
			title=L"\mbox{Total games tanked by cutoff and tanking probability}",
			xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
			yticks=(Array(miny:inc:maxy),[@sprintf("\$%.0f\$", i) for i in miny:inc:maxy]),
			legend=:bottomright,
			#legend=:best,
			legendfont=6,
			legendtitle=L"\mbox{Draft ranking cutoff}",
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
	savefig(fig, string(results_dir,"/avg_games_tanked",".pdf"));

	## Plot avg_already_tank
	print("Plotting avg_already_tank: average number of tanking teams\n")
	miny = Int(ceil(findmin(avg_already_tank)[1]))
	maxy = Int(floor(findmax(avg_already_tank)[1]))
	fig = Plots.plot(show=false,
			xlab=L"\mbox{Probability of tanking once eliminated}", 
			ylab=L"\mbox{Average number of tanking teams}",
			title=L"\mbox{Number of tanking teams by cutoff and tanking probability}",
			xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:10:100]),
			yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
			legend=:topleft,
			legendfont=6,
			legendtitle=L"\mbox{Draft ranking cutoff}",
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
	savefig(fig, string(results_dir,"/avg_already_tank",".pdf"));

	## Plot avg_eliminated
	print("Plotting avg_eliminated: average number of eliminated teams by every game of the season\n")
	miny = Int(ceil(findmin(avg_eliminated)[1]))
	maxy = Int(floor(findmax(avg_eliminated)[1]))
	inc = num_games / 5
	fig = Plots.plot(show=false,
			xlab=L"\mbox{Percent through season}", 
			ylab=L"\mbox{Number teams eliminated}",
			title=L"\mbox{Number teams eliminated over time}",
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
	savefig(fig, string(results_dir,"/avg_eliminated",".pdf"));
end # main
