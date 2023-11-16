#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
feature_type_list=("efficacy" "safety" "combinedEfficacySafety" "kegg" "smpdbDrugMet" "smpdbDrugAct" "misc")



# # Prepare drug combinations
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=DrugComb_$disease --output=Trash/sampo_out/DrugComb_$disease.out --export=disease=$disease sampo/1_sampo_Rscript.sh $disease
#   done

  
  
# # Assigning categories for drug combinations
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=DrugComb_cat_effVadv_$disease --output=Trash/sampo_out/DrugComb_cat_effVadv_$disease.out --export=disease=$disease sampo/7_sampo_Rscript.sh $disease
#   done



# Scatter plot of top features
for feature_type in ${feature_type_list[@]}
  do
    sbatch --job-name=NES_scatter_$feature_type --output=Trash/sampo_out/NES_scatter_$feature_type.out --export=feature_type=$feature_type sampo/9_sampo_Rscript.sh $feature_type
  done
  
  
  
# # Split data for ML
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=ML_dataSplit_$disease --output=Trash/sampo_out/ML_dataSplit_$disease.out --export=disease=$disease sampo/8_sampo_Rscript.sh $disease
#   done


