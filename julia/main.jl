# Needs Plots, DelimitedFiles, GR
using Plots
using DelimitedFiles
using LaTeXStrings
include("tanking.jl")

function main(do_simulation = true, num_repeats = 100000, results_dir = "../results")
	## Variables that need to be set
	num_teams = 30 # number of teams
	num_rounds = 3 # a round consists of each team playing each other team
	num_steps = 20 # discretization of [0,1] for tanking probability
	gamma = 0.75 # probability a better-ranked team wins over a worse-ranked team

	## When the (draft) ranking will be set as fraction of number games
	set_ranking = [4/8; 5/8; 6/8; 7/8; 1]
	num_rankings = length(set_ranking)
	#color_for_cutoff_point = ["c" "b" "m" "r" "k"]
	## end variables that need to be set

	## For output
	avg_kend = 0 # holds KT distance for each tanking probability and cutoff for draft ranking
	avg_already_tank = 0

	if do_simulation
		## Do simulation
		avg_kend, avg_already_tank = simulate(num_teams, num_rounds, num_repeats, num_steps, gamma, set_ranking)
		writedlm(string(results_dir, "/avg_kend.csv"), avg_kend, ',')
		writedlm(string(results_dir, "/avg_already_tank.csv"), avg_already_tank, ',')
	else
		avg_kend = readdlm(string(results_dir, "/avg_kend.csv"), ',')
		avg_already_tank = readdlm(string(results_dir, "/avg_already_tank.csv"), ',')
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

	miny = Int(ceil(findmin(avg_kend)[1]))
	maxy = Int(floor(findmax(avg_kend)[1]))

	fig = Plots.plot(show=false,
			xlab=L"\mbox{Probability of tanking}", 
			ylab=L"\mbox{Kendall tau distance for non-playoff teams}",
			title=L"\mbox{Kendall tau distance}",
			#title=L"\mbox{Kendall tau distance by cutoff and tanking probability}",
			xticks=(Array(0:0.1:1),["\$$i\$" for i in 0:0.1:1]),
			yticks=(Array(miny:maxy),["\$$i\$" for i in miny:maxy]),
			legend=:bottomright,
			legendfont=8,
			titlefont=12,
			tickfont=8,
			grid=false,
			display_type=:inline);
	for r = 1:num_rankings
		cutoff_game = set_ranking[r]
		plot!(0:(1/num_steps):1, avg_kend[:,r], 
				label=latexstring("$cutoff_game", "\\mbox{ of season}"));
	end
	savefig(fig, string(results_dir,"/plot",".pdf"));
	#savefig(fig, "plot.pdf");
	#save(string(results_dir,"/plot",".pdf"), fig);
end # main
