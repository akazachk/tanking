using PackageCompiler

using Combinatorics
using DelimitedFiles
using Distributions
using JuMP
using LaTeXStrings
using Plots
using Printf
using PyCall
using Random
using StatsPlots

pkg = [
       :Combinatorics,
       :DelimitedFiles,
       :Distributions,
       :JuMP,
       :LaTeXStrings,
       :Plots,
       :Printf,
       :PyCall,
       :Random,
       :StatsPlots
      ]

# Gurobi is optional (only needed for math_elim_mode with abs value >= 2)
try
  @eval using Gurobi
  push!(pkg, :Gurobi)
catch
  @warn "Gurobi could not be loaded; building the system image without it"
end

@info "Building system image..."
PackageCompiler.create_sysimage(
  pkg,
  project=".",
  sysimage_path="build/JuliaTanking.so",
  precompile_statements_file="build/precompile.jl"
)
