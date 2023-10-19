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
#SBATCH --ntasks 1   # Number of task
#SBATCH --time 01:00:00   # Runtime
#SBATCH --mem=10000   # Reserve 10 GB RAM for the job
#SBATCH --partition serial   # Partition to submit
#SBATCH --mail-user arindam.ghosh@uef.fi      # this is the email you wish to be notified at
#SBATCH --mail-type ALL   # ALL will alert you of job beginning, completion, failure etc





/opt/R/4.2.2/bin/R --version
/opt/R/4.2.2/bin/Rscript --version
grep -c ^processor /proc/cpuinfo

echo "Prepare drug combinations for " $1 " ---------------------"

/opt/R/4.2.2/bin/Rscript Scripts/Drug_combinations/Extract_drugCombs_and_targets.R --disease $1