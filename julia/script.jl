#!/usr/bin/env julia

import Pkg
Pkg.activate(".")
#Pkg.init()

include("main.jl"); @time main_simulate(do_simulation=true,num_replications=100,do_plotting=false,num_steps=30,math_elim_mode=-2, gamma=0.71425, results_dir="../results/tmp", selected_steps=[1])
