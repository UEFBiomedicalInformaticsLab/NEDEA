#!/bin/bash

disease_list=("LungCancer" "BreastCancer" "ProstateCancer" "OvaryCancer" "KidneyCancer" "SkinCancer")
data_balance_method=("none")
metric_list=("BalancedAccuracy" "F1" "PRAUC")

  


# ML RWR FGSEA
for disease in ${disease_list[@]}
do
  for balance in ${data_balance_method[@]}
  do
    sbatch --job-name=ML_RWR_FGSEA_$disease\_$balance --output=ML_RWR_FGSEA_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/11_sampo_Rscript.sh $disease $balance
  done
done
  
  
  
# ML BBSI Drug-ADR
for balance in ${data_balance_method[@]}
do
  for disease in ${disease_list[@]}
  do
    sbatch --job-name=ML_BBSI_DrgAdr_$disease\_$balance --output=ML_BBSI_DrgAdr_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/12_sampo_Rscript.sh $disease $balance
  done
done



# ML BBSI Drug-Disease
for balance in ${data_balance_method[@]}
do
  for disease in ${disease_list[@]}
  do
    sbatch --job-name=ML_BBSI_DrgDis_$disease\_$balance --output=ML_BBSI_DrgDis_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/13_sampo_Rscript.sh $disease $balance
  done
done



# ML BBSI Drug-Drug
for disease in ${disease_list[@]}
do
  for balance in ${data_balance_method[@]}
  do
    sbatch --job-name=ML_BBSI_DrgDrg_$disease\_$balance --output=ML_BBSI_DrgDrg_$disease\_$balance.out --export=disease=$disease,balance=$balance sampo/14_sampo_Rscript.sh $disease $balance
  done
done




# # Plot train vs test accuracy
# for disease in ${disease_list[@]}
# do
#   for metric in ${metric_list[@]}
#   do
#     sbatch --job-name=plot_trainTestAcc_$disease\_$metric --output=plot_trainTestAcc_$disease\_$metric.out --export=disease=$disease,metric=$metric sampo/19_sampo_Rscript.sh $disease $metric
#   done
# done



# # Plot feature importance
# for disease in ${disease_list[@]}
# do
#   for balance in ${data_balance_method[@]}
#   do
#     sbatch --job-name=plotVarImp_$disease\_$balance --output=plotVarImp_$disease\_$balance.out --export=disease=$disease,balance=$balance, sampo/16_sampo_Rscript.sh $disease $balance
#   done
# done