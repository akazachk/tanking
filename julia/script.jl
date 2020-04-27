#!/usr/bin/env julia
# Run with julia --project

include("Tanking.jl")
#using Tanking

"""
Function main to run simulation with selected steps
"""
function main(selected_steps)
  @time Tanking.main_simulate(do_simulation=true,num_replications=100,do_plotting=false,num_steps=30,math_elim_mode=-2, gamma=0.71425, results_dir="../results/tmp", selected_steps=selected_steps)
end # main

if abspath(PROGRAM_FILE) == @__FILE__
  main(ARGS[1])
end
