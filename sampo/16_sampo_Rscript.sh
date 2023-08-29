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
#SBATCH --ntasks 5   # Number of task
#SBATCH --time 1-00:00:00   # Runtime
#SBATCH --mem=150000   # Reserve 150 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type ALL   # ALL will alert you of job beginning, completion, failure etc





~/miniconda3/envs/interactome/bin/R --version
~/miniconda3/envs/interactome/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "Plot feature importance: " $1 " - " $2, " ---------------------"

~/miniconda3/envs/interactome/bin/Rscript Scripts/Summary_plots/plot_DrugCombs_featureImportance.R --disease $1 --data_balance_method $2