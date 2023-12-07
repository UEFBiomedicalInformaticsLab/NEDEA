#!/bin/bash
#
#
# Name: sampo_Rscript.sh
# Desc.: Run Rscripts in SAMPO
# Updated: 19 October 2023
#
#
#SBATCH --exclusive
#SBATCH --distribution=cyclic
#SBATCH --ntasks 40   # Number of task
#SBATCH --time 01:00:00   # Runtime
#SBATCH --mem=10000   # Reserve 10 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type ALL   # ALL will alert you of job beginning, completion, failure etc





/opt/R/4.2.2/bin/R --version
/opt/R/4.2.2/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "Extend efficacy safety library for " $1 " ---------------------"
/opt/R/4.2.2/bin/Rscript Extend_Eff_Saf/Build_EfficacySafety_lib.R --disease $1 --nproc 20




# echo "Running FGSEA (extended-efficacy-safety) for " $1 "- " $2 " ---------------------"
# /opt/R/4.2.2/bin/Rscript Extend_Eff_Saf/Run_FGSEA.R --disease $1 --drug_target_type $2