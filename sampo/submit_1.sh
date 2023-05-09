#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")



# # Prepare drug combinations
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=DrugComb_$disease --output=DrugComb_$disease.out --export=disease=$disease sampo/1_sampo_Rscript.sh $disease
#   done



# # Running RWR for drug combinations
# for disease in ${disease_list[@]}
#  do
#    sbatch --job-name=RWR_$disease --output=RWR_$disease.out --export=disease=$disease sampo/2_sampo_Rscript.sh $disease
#  done



# # RWR gene ranking check
# for disease in ${disease_list[@]}
#   do
#      sbatch --job-name=RWR_geneRank_$disease --output=RWR_geneRank_$disease.out --export=disease=$disease sampo/3_sampo_Rscript.sh $disease
#    done


  
# # FGSEA common library
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=FGSEA_Common_$disease --output=FGSEA_Common_$disease.out --export=disease=$disease sampo/4_sampo_Rscript.sh $disease
#   done



# # FGSEA efficacy-safety library
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=FGSEA_EffSaf_$disease --output=FGSEA_EffSaf_$disease.out --export=disease=$disease sampo/5_sampo_Rscript.sh $disease
#   done
  
  
  
# # BBSI_DrgAdr
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=BBSI_DrgAdr_$disease --output=BBSI_DrgAdr_$disease.out --export=disease=$disease sampo/6_sampo_Rscript.sh $disease
#   done
  
  
  
# # BBSI_DrgDis
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=BBSI_DrgDis_$disease --output=BBSI_DrgDis_$disease.out --export=disease=$disease sampo/7_sampo_Rscript.sh $disease
#   done
  
  
  
# # BBSI_DrgDrg
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=BBSI_DrgDrg_$disease --output=BBSI_DrgDrg_$disease.out --export=disease=$disease sampo/8_sampo_Rscript.sh $disease
#   done
  
  
    
# # Steiner tree (ADR-Disease-Drug)
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=SteinerTree_AdrDisDrg_$disease --output=SteinerTree_AdrDisDrg_$disease.out --export=disease=$disease sampo/9_sampo_Rscript.sh $disease
#   done
  
  
  
# # ML data split
# for disease in ${disease_list[@]}
#   do
#     sbatch --job-name=ML_dataSplit_$disease --output=ML_dataSplit_$disease.out --export=disease=$disease sampo/10_sampo_Rscript.sh $disease
#   done



# Model accuracy
for disease in ${disease_list[@]}
  do
    sbatch --job-name=plot_ModelAcc_$disease --output=plot_ModelAcc_$disease.out --export=disease=$disease sampo/17_sampo_Rscript.sh $disease
  done
  
  
  
# Model hyperparameters
for disease in ${disease_list[@]}
  do
    sbatch --job-name=plot_hyperParam_$disease --output=plot_hyperParam_$disease.out --export=disease=$disease sampo/18_sampo_Rscript.sh $disease
  done