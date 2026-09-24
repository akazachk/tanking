# On Tanking and Competitive Balance
##### By Aleksandr M. Kazachkov and Shai Vardi
##### February 2020

This project contains the code for a simulator of an NBA season aimed at understanding tanking behavior.

### Requirements
For general requirements, check the "compat" section under [`Project.toml`](Project.toml). The checked-in `Manifest.toml` is in the pre-1.7 format, so use Julia 1.6 (tested with 1.6.7). If `PyPlot` is available, the figures in the paper can be plotted using `do_plotting=true` in the commands below; the pinned `PyPlot` needs matplotlib < 3.9 (e.g., `pip install "matplotlib<3.9"` and set `PYTHON=python3` before instantiating).

On a Mac, for plotting, one needs to install XQuartz.

Gurobi (9.0 or 9.1, for the pinned Gurobi.jl 0.9) is only needed to solve the MIPs for mathematical elimination, i.e., when `abs(math_elim_mode) >= 2` (including the default `math_elim_mode=-2` of `main_simulate`); it is loaded the first time it is needed. Everything else (parsing NBA data, model validation, noisy rankings, and simulations with `math_elim_mode` in -1, 0, 1) runs without it; in that case, `Pkg.instantiate()` reports that Gurobi failed to build/precompile, which can be ignored. To install Gurobi.jl, use `Pkg.build("Gurobi")` in a shell in which the `GUROBI_HOME` variable is defined or `Gurobi` can be found on the `PATH`. E.g., on Mac, `GUROBI_HOME` is set to `/Library/<gurobiversion>/mac64`.

Before running the code, you should [instatiate the environment](https://pkgdocs.julialang.org/v1/environments/). It is further strongly recommended to create a sysimage, the steps for which should be automatically performed if you type `make` from the main project directory on a Linux or Mac.

### Running the code
To run a simulation, parse NBA data, and reproduce data regarding noisiness of the reverse order ranking, start `julia` from the main project directory, type `]` to enter `pkg` mode, and type `activate ./` to activate the tanking environment. Type `instantiate` to get the required packages. Afterwards, pressing `backspace` will return you to the normal prompt. You can also run `julia --project` to avoid the `activate` step above.

The code can be run with the following commands:
				
		using Tanking
		Tanking.main_simulate(do_simulation=1, num_replications=100000, do_plotting=false, mode=Tanking.STRICT, math_elim_mode=-2, gamma=0.71425) 
		Tanking.main_parse(do_plotting=false, mode=Tanking.STRICT)
		Tanking.rankings_are_noisy(do_simulation=true, num_replications=100000, do_plotting=false, mode=Tanking.STRICT)
				
All of the experiments (model validation, simulation, parsing NBA data, noisiness of rankings) can also be rerun with a single script, from the main project directory:

		julia --project=. scripts/run_experiments.jl --results-dir=results/rerun
		julia --project=. scripts/run_experiments.jl --replications=100 --results-dir=results/tmp   # quick test
		julia --project=. scripts/run_experiments.jl --math-elim-mode=0 --results-dir=results/rerun   # without Gurobi

Run `head -35 scripts/run_experiments.jl` to see all options (e.g., `--gamma=auto` to use the value of gamma that best fits the NBA data, or `--steps` / `--aggregate` to split the simulation across jobs).

### NBA data
The directory [`data`](data) contains the results of every regular-season game (from [basketball-reference.com](https://www.basketball-reference.com)) in `data/gamesYYZZ.csv` for the seasons 2004-05 through 2025-26, except 2011-12 (lockout) and 2019-20 and 2020-21 (COVID-19), in which teams did not play 82 games. The list of seasons that is used is `Tanking.nba_seasons`; pass `seasons=Tanking.nba_seasons_2004_2019` to `main_parse` or `BT_MLE` to use only the seasons in the original paper. The file `data/winpct.csv` contains the win percentage of the team in each final position (rows) for every season (columns).

The files for 2021-22 through 2025-26 were downloaded in September 2026; every team's record in them matches the basketball-reference standings. To (re)download seasons and regenerate `data/winpct.csv` (a season is named by the year in which it ends):

		python3 scripts/fetch_bbref_games.py 2022 2023 2024 2025 2026

Play-in games, playoff games, and the NBA Cup championship game (which does not count in the standings) are excluded. Note that `data/winpct.csv` uses each team's actual number of games, so the 2012-13 column differs slightly from the original file for Boston and Indiana (81 games, after their canceled game).


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
