#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
<<<<<<< HEAD
data_balance_method=("none")
models=("rf" "glmnet" "svmRadial" "nb")
=======
data_balance_method=("none" "SMOTE" "upSample" "downSample")
models=("rf" "glmnet" "svmRadial" "knn" "nb")
>>>>>>> dd66bdd55a3da78129090252ed59959b311d68ad
proximity_list=("separation" "closest" "shortest" "centre" "kernel")

 
 
  
  
<<<<<<< HEAD
# # Plot pan cancer accuracy
# sbatch --job-name=plotPanCancer --output=plotPanCancer.out sampo/20_sampo_Rscript.sh 


# ML BBSI Drug-Disease-ADR 
for balance in ${data_balance_method[@]}
do
  for proximity in ${proximity_list[@]}
  do
    for disease in ${disease_list[@]}
    do
    sbatch --job-name=ML_BBSI_DrgDisAdr_$disease\_$proximity\_$balance --output=ML_BBSI_DrgDisAdr_$disease\_$proximity\_$balance.out --export=disease=$disease,balance=$balance,proximity=$proximity, sampo/23_sampo_Rscript.sh $disease $balance $proximity
    done
  done
done
=======
# Plot pan cancer accuracy
sbatch --job-name=plotPanCancer --output=plotPanCancer.out sampo/20_sampo_Rscript.sh 


# # ML BBSI Drug-Disease-ADR 
# for balance in ${data_balance_method[@]}
# do
#   for proximity in ${proximity_list[@]}
#   do
#     for disease in ${disease_list[@]}
#     do
#     sbatch --job-name=ML_BBSI_DrgDisAdr_$disease\_$proximity\_$balance --output=ML_BBSI_DrgDisAdr_$disease\_$proximity\_$balance.out --export=disease=$disease,balance=$balance,proximity=$proximity, sampo/23_sampo_Rscript.sh $disease $balance $proximity
#     done
#   done
# done
>>>>>>> dd66bdd55a3da78129090252ed59959b311d68ad


