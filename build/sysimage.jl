using PackageCompiler

pkg = [
       :Combinatorics,
       :DelimitedFiles,
       :Distributions,
       :Gurobi,
       :JuMP,
       :LaTeXStrings,
       :MathOptInterface,
       :Printf,
       :PyCall,
       :PyPlot,
       :Random,
       :Tanking
      ]

@info "Building system image..."
PackageCompiler.create_sysimage(
  pkg,
  project=".",
  sysimage_path="build/JuliaTanking.so",
  precompile_statements_file="build/precompile.jl"
)
