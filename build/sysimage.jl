using PackageCompiler

using Combinatorics
using DelimitedFiles
using Distributions
using Gurobi
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
       :Gurobi,
       :JuMP,
       :LaTeXStrings,
       :Plots,
       :Printf,
       :PyCall,
       :Random,
       :StatsPlots
      ]

@info "Building system image..."
PackageCompiler.create_sysimage(
  pkg,
  project=".",
  sysimage_path="build/JuliaTanking.so",
  precompile_statements_file="build/precompile.jl"
)
