# On Tanking and Competitive Balance
##### By Aleksandr M. Kazachkov and Shai Vardi
##### February 2020

This project contains the code for a simulator of an NBA season aimed at understanding tanking behavior.

### Requirements
For general requirements, check the "compat" section under [`Project.toml`](Project.toml). If `PyPlot` is available, the figures in the paper can be plotted using `do_plotting=true` in the commands below.

On a Mac, for plotting, one needs to install XQuartz.

Gurobi is needed. To install, you need to use `Pkg.build("Gurobi")` in a shell in which the `GUROBI_HOME` variable is defined or `Gurobi` can be found on the `PATH`. E.g., on Mac, `GUROBI_HOME` is set to `/Library/<gurobiversion>/mac64`.

Before running the code, you should [instatiate the environment](https://pkgdocs.julialang.org/v1/environments/). It is further strongly recommended to create a sysimage, the steps for which should be automatically performed if you type `make` from the main project directory on a Linux or Mac.

### Running the code
To run a simulation, parse NBA data, and reproduce data regarding noisiness of the reverse order ranking, first change directories to the `julia` subdirectory, then start `julia`, type `]` to enter `pkg` mode, and type `activate ./` to activate the tanking environment. Type `instantiate` to get the required packages. Afterwards, pressing `backspace` will return you to the normal prompt. You can also run `julia --project` to avoid the `activate` step above.

The code can be run with the following commands:
				
		using Tanking
		Tanking.main_simulate(do_simulation=1, num_replications=100000, do_plotting=false, mode=Tanking.STRICT, math_elim_mode=-2, gamma=0.71425) 
		Tanking.main_parse(do_plotting=false, mode=Tanking.STRICT)
		Tanking.rankings_are_noisy(do_simulation=true, num_replications=100000, do_plotting=false, mode=Tanking.STRICT)
				

### Options
1. Option `mode` repesents the base true ranking.
				
		mode = STRICT: true ranking is strict. 1 > 2 > ... > 30
		mode = TIES: true ranking has ties. [1,5] > [6,10] > ... > [26,30]
		mode = BT_DISTR: each team gets a strength score based on random distribution (either Beta(1,1) = uniform, or Beta(2,5)), and game winners are determined by (Zermelo-)Bradley-Terry model
		mode = BT_EXPONENTIAL: same as mode = 3, except winners are determined with exponential version of the Bradley-Terry model
		mode = BT_ESTIMATED: each team gets a strength score (calibrated using an MLE on NBA data) and game winners are determined by (Zermelo-)Bradley-Terry model
				
2. Number of teams can be changed in the code (`num_teams`).
3. There are other plotting mechanisms implemented, but not all have been tested thoroughly.
4. `math_elim_mode`: how an eliminated team is identified, and how this will be used to determine tanking

    0: use effective elimination
    1: use mathematical elimination, but calculated by heuristics only
    2: use math elim, binary MIP, team-wise formulation
    3: use math elim, general integer MIP, team-wise formulation
    4: use math elim, binary MIP, cutoff formulation
    5: use math elim, general integer MIP, cutoff formulation
    <0: use effective elimination for tanking, but calculate mathematical elimination

### Assumptions
1. No simulataneous games
2. No conference / division play
3. True ranking is static
4. No home/away advantage
5. Ties broken between two teams with the same win percentage is by head-to-head record, and afterwards uniformly at random


### To create a sysimage

1. Create pre-compilation statements (running from project directory)

        mkdir -p results/tmp
        julia --trace-compile="precompile.jl" --project="Tanking" scripts/test_script.jl 1

2. Create the sysimage (from the Tanking directory)

        julia> using PackageCompiler
        julia> PackageCompiler.create_sysimage([:Combinatorics, :DelimitedFiles, :Distributions, :Gurobi, :JuMP, :LaTeXStrings, :MathOptFormat, :Plots, :Printf, :PyCall, :Random, :StatsPlots, :Tanking], project=".", sysimage_path="JuliaTanking.so", precompile_statements_file="precompile.jl")

3. Run the script

        julia --sysimage=Tanking/JuliaTanking.so --project="Tanking" scripts/run_script.jl
