#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 96:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 100   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "s6.FC"   # job name
#SBATCH --mem=300G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=hugj2006@iastate.edu   # email address
#SBATCH --partition=whatever -- use sinfo to get queue names


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
module load r/3.5.0-py2-ufvuwmm
R CMD BATCH 6b.FC.r
