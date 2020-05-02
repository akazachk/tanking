#!/bin/bash
#SBATCH --array=1
#SBATCH --time=36:00:00
#SBATCH --account=def-alodi
#SBATCH --mem-per-cpu=32G
#SBATCH --cpus-per-task=1
julia --sysimage=/home/akazachk/projects/def-alodi/akazachk/tanking/Tanking/JuliaTanking.so --project=/home/akazachk/projects/def-alodi/akazachk/tanking/Tanking/ /home/akazachk/projects/def-alodi/akazachk/tanking/scripts/run_script.jl $SLURM_ARRAY_TASK_ID /home/akazachk/projects/def-alodi/akazachk/tanking/results/2020-04-30
