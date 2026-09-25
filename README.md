# On Tanking and Competitive Balance
##### By Aleksandr M. Kazachkov and Shai Vardi
##### February 2020

This project contains the code for a simulator of an NBA season aimed at understanding tanking behavior.

### Requirements
For general requirements, check the "compat" section under [`Project.toml`](Project.toml). The checked-in `Manifest.toml` was resolved with **Julia 1.13** (JuMP 1.x, MathOptInterface 1.x, Gurobi.jl 1.9 with Gurobi 13); use Julia 1.13 or later (with [juliaup](https://github.com/JuliaLang/juliaup): `juliaup add 1.13`, then `julia +1.13`, or `JULIA="julia +1.13" scripts/run_all.sh ...`). Gurobi.jl uses the Gurobi 13 library from the `Gurobi_jll` package, so no separate Gurobi installation is needed, only a license (e.g., a `gurobi.lic` file found by Gurobi, or `GRB_LICENSE_FILE` pointing to it); to use a local Gurobi installation instead, set `GUROBI_JL_USE_GUROBI_JLL=false` and `GUROBI_HOME` and run `Pkg.build("Gurobi")`. The code as of the 2020 experiments used Julia 1.6 and JuMP 0.21 (see the git history). If `PyPlot` is available, the figures in the paper can be plotted using `do_plotting=true` in the commands below; set `PYTHON=python3` (a Python with matplotlib) before instantiating so that PyCall uses it.

On a Mac, for plotting, one needs to install XQuartz.

Gurobi (9.0 or 9.1, for the pinned Gurobi.jl 0.9) is only needed to solve the MIPs for mathematical elimination, i.e., when `abs(math_elim_mode) >= 2` (including the default `math_elim_mode=-2` of `main_simulate`); it is loaded the first time it is needed. Everything else (parsing NBA data, model validation, noisy rankings, and simulations with `math_elim_mode` in -1, 0, 1) runs without it; in that case, `Pkg.instantiate()` reports that Gurobi failed to build/precompile, which can be ignored. To install Gurobi.jl, use `Pkg.build("Gurobi")` in a shell in which the `GUROBI_HOME` variable is defined or `Gurobi` can be found on the `PATH`. E.g., on Mac, `GUROBI_HOME` is set to `/Library/<gurobiversion>/mac64`.

Before running the code, you should [instantiate the environment](https://pkgdocs.julialang.org/v1/environments/). It is further strongly recommended to create a sysimage, the steps for which should be automatically performed if you type `make` from the main project directory on a Linux or Mac.

### Running the code
To run a simulation, parse NBA data, and reproduce data regarding noisiness of the reverse order ranking, start `julia` from the main project directory, type `]` to enter `pkg` mode, and type `activate ./` to activate the tanking environment. Type `instantiate` to get the required packages. Afterwards, pressing `backspace` will return you to the normal prompt. You can also run `julia --project` to avoid the `activate` step above.

The code can be run with the following commands:
				
		using Tanking
		loss, gammas, gamma = Tanking.model_validation(num_replications=100000, do_plotting=false)   # gamma chosen by the minimax rule
		Tanking.main_simulate(do_simulation=1, num_replications=100000, do_plotting=false, mode=Tanking.STRICT, math_elim_mode=-2, gamma=gamma) 
		Tanking.main_parse(do_plotting=false, mode=Tanking.STRICT)
		Tanking.rankings_are_noisy(do_simulation=true, num_replications=100000, do_plotting=false, mode=Tanking.STRICT)
				
### Recreating all experiments and plots
`scripts/run_all.sh` reruns everything (model validation, the simulation, NBA data, noisy rankings, and the sensitivity run below), with plots, from the main project directory:

		test/test_all.sh                              # quick test first (20 replications, all seasons, needs Gurobi)
		scripts/run_all.sh -o results/final -j 16     # full run: 100K replications, all seasons, 16 parallel jobs

`test/test_all.sh` runs `scripts/run_all.sh` with few replications (including splitting the simulation over parallel jobs and aggregating it) and then checks every output with `test/check_results.jl` (file sizes, finite values, min <= avg <= max, that MIPs were solved when Gurobi is needed, that the requested NBA seasons were used, and the sensitivity checks); it ends with `PASSED` or `FAILED`. Both scripts take the same options (run with `-h`), e.g., `-s 2004-2019` for the 14 seasons used in the 2020 experiments, `-g <value>` to fix gamma instead of using the value chosen by model validation (e.g., `-g 0.71425`, the value used in the 2020 experiments) (the default, `-g auto`: the value in the grid whose largest error over 0, 15, and 30 selfish teams is smallest; it is saved in `gamma.txt`, and a warning is shown if the grid is too coarse around it), `-m 0` to run without Gurobi, `-N` for no plots. Plots need LaTeX (with `dvipng`). Logs of every step go to `<results dir>/logs`, and the settings to `<results dir>/settings.txt`.

`test/test_mip.jl` tests the mathematical-elimination MIPs without a Gurobi license, using the open-source solver HiGHS (`julia --project=test test/test_mip.jl`, from the main project directory; the `test` environment adds HiGHS): every elimination check is solved with both exact formulations (binary and general integer), which must agree, and the stored best schedules are checked for consistency.

Model validation compares the simulated win percentage by rank with the NBA data for 0, 15 (on average: each team is selfish with probability 1/2), and 30 selfish teams (plots `win_pct_{0,15,30}selfish`, and the model-error plot `model_loss`). `settings.txt` in the results directory records the start and finish times, the duration, the git commit, and any uncommitted changes to the code.

With `-j` > 1, model validation (one job per model) and the simulation (one job per step) run in parallel; the noisy-rankings experiment is a single process, so time it with a test run first (e.g., `-n 100`) and, if needed, run it separately with fewer replications (`-e noisy -n ...`).

The individual experiments can also be run with `scripts/run_experiments.jl`:

		julia --project=. scripts/run_experiments.jl --results-dir=results/rerun
		julia --project=. scripts/run_experiments.jl --math-elim-mode=0 --results-dir=results/rerun   # without Gurobi

Run `head -40 scripts/run_experiments.jl` to see all options (e.g., `--seasons`, `--gamma=auto`, or `--steps` / `--aggregate` to split the simulation across jobs).

### Sensitivity to tanking after the breakpoint
By default, a selfish team tanks in every game after it is eliminated, including games after the breakpoint, so that one simulated season serves all breakpoints. `scripts/run_sensitivity.jl` measures the effect of this on the Kendall tau results: on seasons that are identical up to the breakpoint (common random numbers), it compares the default behavior with one in which no team tanks after the breakpoint (`simulate(...; stop_tanking_after_breakpoint=true, seed_per_replication=seed)`), using effective elimination (`math_elim_mode = 0`, which makes the same tanking decisions as the default `-2`, without Gurobi):

		julia --project=. -t 4 scripts/run_sensitivity.jl 10000 results/sens_test "[1,9,16,24,31]"

Step `s` corresponds to `s-1` selfish teams (default: all 31 steps). The script prints two checks that must be exactly 0 (games tanked up to each breakpoint, and the Kendall tau at the end of the season, are the same in both behaviors), and writes `kend_keep.csv`, `kend_stop.csv`, `kend_diff.csv` (stop minus keep), standard errors `se_*.csv` (for the difference, computed from the paired differences), `games_tanked_*.csv`, and `breakpoints.csv`. Because each replication is seeded separately, the baseline does not reproduce the 2020 experiment results (in `results/2020-*`) exactly, only statistically.

### NBA data
The directory [`data`](data) contains the results of every regular-season game (from [basketball-reference.com](https://www.basketball-reference.com)) in `data/gamesYYZZ.csv` for the seasons 2004-05 through 2025-26, except 2011-12 (lockout) and 2019-20 and 2020-21 (COVID-19), in which teams did not play 82 games. The list of seasons that is used is `Tanking.nba_seasons`; pass `seasons=Tanking.nba_seasons_2004_2019` to `main_parse` or `BT_MLE` to use only the seasons used in the 2020 experiments. The file `data/winpct.csv` contains the win percentage of the team in each final position (rows) for every season (columns).

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
