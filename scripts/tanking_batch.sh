#!/bin/bash
#SBATCH --array=1-31
#SBATCH --time=00:01:00
#SBATCH --account=def-alodi
#SBATCH --mem-per-cpu=16G
#SBATCH --cpus-per-task=1
julia --sysimage=/home/akazachk/projects/def-alodi/akazachk/tanking/julia/JuliaTanking.so --project=/home/akazachk/projects/def-alodi/akazachk/tanking/julia/ /home/akazachk/projects/def-alodi/akazachk/tanking/scripts/script.jl $SLURM_ARRAY_TASK_ID /home/akazachk/projects/def-alodi/akazachk/tanking/results/tmp
