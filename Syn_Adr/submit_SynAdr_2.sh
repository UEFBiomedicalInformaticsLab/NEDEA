#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
drug_target_type_list=("known" "PS" "SIGNOR" "NPA" "RI" "KEGG" "all")


# Execute FGSEA (allADR) on drug combinations
for disease in ${disease_list[@]}
do
  for drug_target_type in ${drug_target_type_list[@]}
  do
    sbatch --job-name=synAdr_FGSEA_allAdr_$disease\_$drug_target_type --output=Trash/sampo_out/synAdr_FGSEA_allAdr_$disease\_$drug_target_type.out --export=disease=$disease,drug_target_type=$drug_target_type Syn_Adr/sampo_SynAdr_1.sh $disease $drug_target_type
  done
done

