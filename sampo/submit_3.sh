#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
data_balance_method=("none" "SMOTE" "upSample" "downSample")
models=("rf" "glmnet" "svmRadial" "knn" "nb")


 
# Plot feature importance  
for disease in ${disease_list[@]}
do
  for balance in ${data_balance_method[@]}
  do
    for model in ${models[@]}
    do
    sbatch --job-name=plotVarImp_$disease\_$model\_$balance --output=plotVarImp_$disease\_$model\_$balance.out --export=disease=$disease,balance=$balance,model=$model, sampo/16_sampo_Rscript.sh $disease $balance $model
    done
  done
done
  
  
  
# Plot pan cancer accuracy
sbatch --job-name=plotPanCancer --output=plotPanCancer.out sampo/20_sampo_Rscript.sh 