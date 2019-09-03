# Reducing Tanking Incentives in the NBA
##### By Aleksandr M. Kazachkov and Shai Vardi
##### December 2018

This project contains the code for a simulator of an NBA season aimed at understanding tanking behavior.

### Requirements
Julia 1.x. If `PyPlot` is available, the figures in the paper can be plotted using `do_plotting=true` in the commands below.

### Running the code
To run a simulation, parse NBA data, and reproduce data regarding noisiness of the reverse order ranking, first change directories to the `julia` subdirectory, then start `julia`, type `]` to enter `pkg` mode, and type `activate ./` to activate the tanking environment. Type `instantiate` to get the required packages. Afterwards, pressing `backspace` will return you to the normal prompt. Afterwards, the code can be run with the following commands:
				
		include("main.jl")
		main_simulate(do_simulation=true, num_replications=100000, do_plotting=false, mode=1) 
		main_parse(do_plotting=false, mode=1) 
		rankings_are_noisy(do_simulation=true, num_replications=100000, do_plotting=false, mode=1) 
				

### Options
1. Option `mode` repesents the base true ranking.
				
		mode = 1: true ranking is strict. 1 > 2 > ... > 30
		mode = 2: true ranking has ties. [1,5] > [6,10] > ... > [26,30]
		mode = 3: each team gets a strength score from U[0,1] and game winners are determined by (Zermelo-)Bradley-Terry model
		mode = 4: same as mode = 3, except winners are determined with exponential version of the Bradley-Terry model
				
2. Number of teams can be changed in the code (`num_teams`).
3. There are other plotting mechanisms implemented, but not all have been tested thoroughly.

### Assumptions
1. No simulataneous games
2. No conference / division play
3. True ranking is static
4. No home/away advantage
