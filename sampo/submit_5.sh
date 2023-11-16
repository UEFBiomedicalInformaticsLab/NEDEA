#!/bin/bash


disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer") 
drug_target_type_list=("known" "PS" "SIGNOR" "NPA" "RI" "KEGG" "all")
feature_type_list=("efficacy" "safety" "combinedEfficacySafety" "kegg" "smpdbDrugMet" "smpdbDrugAct" "misc")
drugComb_category_type_list=("SL" "SS" "TE" "ME" "ADR" "EA")
top9varying_list=("TRUE" "FALSE")


# Execute FGSEA (pathways) on drug combinations
for disease in ${disease_list[@]}
do
  for drug_target_type in ${drug_target_type_list[@]}
  do
    for feature_type in ${feature_type_list[@]}
    do
      for drugComb_category_type in ${drugComb_category_type_list[@]}
      do
        for top9varying in ${top9varying_list[@]}
        do
        sbatch --job-name=plot_meanNes_$disease\_$drug_target_type\_$feature_type\_$drugComb_category_type\_$top9varying \
               --output=Trash/sampo_out/plot_meanNes_$disease\_$drug_target_type\_$feature_type\_$drugComb_category_type\_$top9varying.out \
               --export=disease=$disease,drug_target_type=$drug_target_type,$feature_type=$feature_type,$drugComb_category_type=$drugComb_category_type,$top9varying=$top9varying \
               sampo/6_sampo_Rscript.sh $disease $drug_target_type $feature_type $drugComb_category_type $top9varying
        done
      done
    done
  done
done



# sbatch --job-name=plot_meanNes_$disease\_$drug_target_type\_$feature_type\_$drugComb_category_type\_$top9varying --output=Trash/Sampo_out/plot_meanNes_$disease\_$drug_target_type\_$feature_type\_$drugComb_category_type\_$top9varying.out --export=disease=$disease,drug_target_type=$drug_target_type,$feature_type=$feature_type,$drugComb_category_type=$drugComb_category_type,$top9varying=$top9varying sampo/6_sampo_Rscript.sh $disease $drug_target_type $feature_type $drugComb_category_type $top9varying
