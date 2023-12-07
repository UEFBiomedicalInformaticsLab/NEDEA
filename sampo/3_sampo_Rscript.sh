#!/bin/bash
#
#
# Name: sampo_Rscript.sh
# Desc.: Run Rscripts in SAMPO
# Updated: 19 October 2023
#
#
#SBATCH --distribution=cyclic
#SBATCH --ntasks 80   # Number of task
#SBATCH --time 05:00:00   # Runtime
#SBATCH --mem=50000   # Reserve 50 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type ALL   # ALL will alert you of job beginning, completion, failure etc





/opt/R/4.2.2/bin/R --version
/opt/R/4.2.2/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "Running FGSEA (efficacy-safety) for " $1 "- " $2 " ---------------------"

/opt/R/4.2.2/bin/Rscript Scripts/Feature_generation/Execute_FGSEA_EfficacySafety_on_drugCombs.R --disease $1 --drug_target_type $2