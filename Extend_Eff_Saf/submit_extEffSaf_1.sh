#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
feature_type_list=("efficacy" "safety" "combinedEfficacySafety" "kegg" "smpdbDrugMet" "smpdbDrugAct" "misc")
drug_target_type_list=("known" "PS" "SIGNOR" "NPA" "RI" "KEGG" "all")


# Prepare drug combinations
for disease in ${disease_list[@]}
  do
    sbatch --job-name=ExtEnrich_$disease --output=Trash/sampo_out/ExtEnrich_$disease.out --export=disease=$disease Extend_Eff_Saf/sampo_extEffSaf_1.sh $disease
  done


# # Execute FGSEA (extended-efficacy-safety) on drug combinations
# for disease in ${disease_list[@]}
# do
#   for drug_target_type in ${drug_target_type_list[@]}
#   do
#     sbatch --job-name=FGSEA_extEffSaf_$disease\_$drug_target_type --output=Trash/sampo_out/FGSEA_extEffSaf_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=$drug_target_type Extend_Eff_Saf/sampo_extEffSaf_1.sh $disease $drug_target_type
#   done
# done