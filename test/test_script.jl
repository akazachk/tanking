#!/usr/bin/env julia
# Run with julia --project=@Tanking.jl/test

using Tanking

"""
    str2range

Parse string (input should be a single number or range as a string) into UnitRange or into just an Int if it is contains no colon
"""
function str2range(input::AbstractString)
  if findfirst(':', input) != nothing
    y = split(input, ':')
    return UnitRange(parse.(Int,y[1]),parse.(Int,y[2]))
  else
    return parse.(Int,input)
  end
end # str2range


"""
    str2arr

Convert string into array
"""
function str2arr(input::AbstractString)
  # First split the array up if it is given as an array
  y = split(input, x -> (x == '[' || x == ',' || x == ']' || isspace(x)) ? true : false)

  # Remove all empty strings
  y = [y[i] for i in 1:length(y) if y[i] != ""]

  # We now have arr = an array of substrings, each of which should contain an Int _or_ UnitRange
  arr = [str2range(y[i]) for i in 1:length(y)]

  return arr
end # str2arr


"""
    main

Run simulation with selected steps
"""
function main(sel_steps, outdir=nothing)
  if isnothing(outdir) || !isa(outdir,String) || isempty(outdir)
    outdir=ENV["REPOS_DIR"]
    outdir=outdir*"/tanking/results/tmp"
  end
  if !isnothing(sel_steps)
    sel_steps=str2arr(sel_steps)
  end
  @time Tanking.main_simulate(do_simulation=true,num_replications=1,do_plotting=false,num_steps=30,math_elim_mode=-2, gamma=0.71425, results_dir=outdir, selected_steps=sel_steps)
end # main


if abspath(PROGRAM_FILE) == @__FILE__
  if length(ARGS) < 1
    main(nothing)
  elseif length(ARGS) == 1
    main(ARGS[1])
  else
    main(ARGS[1], ARGS[2])
  end
end
