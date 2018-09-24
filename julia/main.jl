using Plots
using DelimitedFiles
include("tanking.jl")

function main(do_simulation = true, results_dir = "../results")
	## Variables that need to be set
	num_repeats = 10000 #10000 #20000

	num_teams = 30 # number of teams
	num_rounds = 3 # a round consists of each team playing each other team
	num_steps = 10 # discretization of [0,1] for tanking probability
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
		writedlm(string(results_dir, "/avg_already_tank.csv"), avg_kend, ',')
	else
		avg_kend = readdlm(string(results_dir, "/avg_kend.csv"), ',')
		avg_already_tank = readdlm(string(results_dir, "/avg_already_tank.csv"), ',')
	end

	## Do plotting
	#pgfplots() # Pkg.add("PGFPlots")
	#pyplot() # Pkg.add("PyPlot") Pkg.add("PyCall") Pkg.add("LaTeXStrings")
	gr() # Pkg.add("GR")

	fig = plot(
			xlab="Probability of tanking", 
			ylab="Kendall tau distance for non-playoff teams", 
			title="Kendall tau distance by cutoff and tanking probability",
			legend=:bottomright,
			grid=false)
	for r = 1:num_rankings
		cutoff_game = set_ranking[r]
		plot!(0:(1/num_steps):1, avg_kend[:,r], label="$cutoff_game of season")
	end
	savefig(fig, string(results_dir,"/plot",".png"))
end # main
