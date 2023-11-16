#!/bin/bash
#
#
# Name: sampo_Rscript.sh
# Desc.: Run Rscripts in SAMPO
# Updated: 19 October 2023
#
#
#SBATCH --ntasks 5   # Number of task
#SBATCH --time 00:30:00   # Runtime
#SBATCH --mem=5000   # Reserve 5 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type FAIL,REQUEUE   # ALL will alert you of job beginning, completion, failure etc





/opt/R/4.2.2/bin/R --version
/opt/R/4.2.2/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "plot the mean NES from FGSEA for " $1 "- " $2 " ---------------------"

/opt/R/4.2.2/bin/Rscript Scripts/Plots/plot_mean_NES_byDrugCategory.R --disease $1 --drug_target_type $2 --feature_type $3 --drugComb_category_type $4 --top9varying $5


