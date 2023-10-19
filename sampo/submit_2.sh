#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
drug_target_type_list=("known" "PS" "SIGNOR" "NPA" "RI" "KEGG" "all")

# data_balance_method=("none")
# metric_list=("BalancedAccuracy" "F1" "PRAUC")
# models=("rf" "svmRadial" "nb" "glmnet")
# feature_types=("Disease2Gene" "WithdrawalAdr2Gene" "CombinedDisAdr2Gene" "keggPath" "SMPDbPath_DrugMet" "SMPDbPath_DrugAction" "miscGeneSet")
  


# # Execute RWR on drug combinations
# for disease in ${disease_list[@]}
# do
#   for drug_target_type in ${drug_target_type_list[@]}
#   do
#     sbatch --job-name=RWR_$disease\_$drug_target_type --output=Trash/Sampo_out/RWR_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=drug_target_type sampo/2_sampo_Rscript.sh $disease $drug_target_type
#   done
# done



# Execute FGSEA (efficacy-safety) on drug combinations
for disease in ${disease_list[@]}
do
  for drug_target_type in ${drug_target_type_list[@]}
  do
    sbatch --job-name=FGSEA_effsaf_$disease\_$drug_target_type --output=Trash/Sampo_out/FGSEA_effsaf_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=drug_target_type sampo/3_sampo_Rscript.sh $disease $drug_target_type
  done
done
  
  
  
# # Execute FGSEA (pathways) on drug combinations
# for disease in ${disease_list[@]}
# do
#   for drug_target_type in ${drug_target_type_list[@]}
#   do
#     sbatch --job-name=FGSEA_path_$disease\_$drug_target_type --output=Trash/Sampo_out/FGSEA_path_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=drug_target_type sampo/4_sampo_Rscript.sh $disease $drug_target_type
#   done
# done  



# # Execute FGSEA (misc) on drug combinations
# for disease in ${disease_list[@]}
# do
#   for drug_target_type in ${drug_target_type_list[@]}
#   do
#     sbatch --job-name=FGSEA_misc_$disease\_$drug_target_type --output=Trash/Sampo_out/FGSEA_misc_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=drug_target_type sampo/5_sampo_Rscript.sh $disease $drug_target_type
#   done
# done