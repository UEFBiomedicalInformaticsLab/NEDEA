#!/bin/bash
#
#
# Name: sampo_Rscript.sh
# Desc.: Run Rscripts in SAMPO
# Updated: 25 April 2023
#
#
#SBATCH --exclusive
#SBATCH --distribution=cyclic
#SBATCH --ntasks 40   # Number of task
#SBATCH --time 03-00:00:00   # Runtime
#SBATCH --mem=50000   # Reserve 50 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type ALL   # ALL will alert you of job beginning, completion, failure etc
#SBATCH --output geneSet_Sim.out	# Standard out goes to this file
#SBATCH --error geneSet_Sim.err	# Standard err goes to this file




~/miniconda3/envs/interactome/bin/R --version
~/miniconda3/envs/interactome/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "GeneSet sim ---------------------"

Rscript Explore_Tox_Dis_sim/GeneSet_Sim_Net.R