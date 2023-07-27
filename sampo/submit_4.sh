#!/bin/bash


# # Final model training
# disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
# data_balance_method=("none")
# models=("rf" "svmRadial" "nb")
# feature_types=("CombinedDisAdr2Gene")
# for disease in ${disease_list[@]}
# do
#  for balance in ${data_balance_method[@]}
#  do
#    for model in ${models[@]}
#    do
#      for featureType in ${feature_types[@]}
#      do
#      sbatch --job-name=ML_FinalModel_$disease\_$balance\_$model\_$featureType --output=ML_FinalModel_$disease\_$balance\_$model\_$featureType.out --export=disease=$disease,balance=$balance,model=$model,featureType=$featureType, sampo/21_sampo_Rscript.sh $disease $data_balance_method $model $featureType
#      done
#    done
#  done
# done



# # External predictions (CDCDB/RX/OTC)
# disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
# data_balance_method=("none")
# models=("rf" "svmRadial" "nb")
# 
# feature_types=("CombinedDisAdr2Gene")
# 
# for disease in ${disease_list[@]}
# do
#   for balance in ${data_balance_method[@]}
#   do
#     for model in ${models[@]}
#     do
#       for featureType in ${feature_types[@]}
#       do
#       sbatch --job-name=ExtPred_CDCDBrxotc_$disease\_$balance\_$model\_$featureType --output=ExtPred_CDCDBrxotc_$disease\_$balance\_$model\_$featureType.out --export=disease=$disease,balance=$balance,model=$model,featureType=$featureType, sampo/24_sampo_Rscript.sh $disease $data_balance_method $model $featureType
#       done
#     done
#   done
# done



# External predictions (CDCDB/AACT)
disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
data_balance_method=("none")
models=("rf" "svmRadial" "nb")
feature_types=("CombinedDisAdr2Gene")

for disease in ${disease_list[@]}
do
  for balance in ${data_balance_method[@]}
  do
    for model in ${models[@]}
    do
      for featureType in ${feature_types[@]}
      do
      sbatch --job-name=ExtPred_CDCDBaact_$disease\_$balance\_$model\_$featureType --output=ExtPred_CDCDBaact_$disease\_$balance\_$model\_$featureType.out --export=disease=$disease,balance=$balance,model=$model,featureType=$featureType, sampo/27_sampo_Rscript.sh $disease $data_balance_method $model $featureType
      done
    done
  done
done



# # External predictions (CDD)
# disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
# data_balance_method=("none")
# models=("rf" "svmRadial" "nb")
# feature_types=("CombinedDisAdr2Gene")

# for disease in ${disease_list[@]}
# do
  # for balance in ${data_balance_method[@]}
  # do
    # for model in ${models[@]}
    # do
      # for featureType in ${feature_types[@]}
      # do
      # sbatch --job-name=ExtPred_CDD_$disease\_$balance\_$model\_$featureType --output=ExtPred_CDD_$disease\_$balance\_$model\_$featureType.out --export=disease=$disease,balance=$balance,model=$model,featureType=$featureType, sampo/26_sampo_Rscript.sh $disease $data_balance_method $model $featureType
      # done
    # done
  # done
# done