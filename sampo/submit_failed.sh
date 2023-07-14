#!/bin/bash


# ML RWR FGSEA
disease_list=("LungCancer" "SkinCancer")
data_balance_method=("none")
for disease in ${disease_list[@]}
do
  for balance in ${data_balance_method[@]}
  do
    sbatch --job-name=ML_RWR_FGSEA_$disease\_$balance --output=ML_RWR_FGSEA_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/11_sampo_Rscript.sh $disease $balance
  done
done



# ML BBSI Drug-Disease
disease_list=("SkinCancer")
data_balance_method=("none")
for balance in ${data_balance_method[@]}
do
  for disease in ${disease_list[@]}
  do
    sbatch --job-name=ML_BBSI_DrgDis_$disease\_$balance --output=ML_BBSI_DrgDis_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/13_sampo_Rscript.sh $disease $balance
  done
done