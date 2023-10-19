#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")



# Prepare drug combinations
for disease in ${disease_list[@]}
  do
    sbatch --job-name=DrugComb_$disease --output=DrugComb_$disease.out --export=disease=$disease sampo/1_sampo_Rscript.sh $disease
  done
