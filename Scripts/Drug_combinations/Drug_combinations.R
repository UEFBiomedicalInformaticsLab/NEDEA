set.seed(5081)


# Script to extract cancer type specific drug combinations



# Read the FIMM drug combinations
FimmDrugComb_drugCombCat <- readRDS("InputFiles/Reference_list/FimmDrugComb_drugCombinations.rds")


# Read the DrugBank drug-drug interactions
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")


# Merging the dataframes based on matching drug ID combinations
FIMM_DrugBank_drugCombs_all <- merge(FimmDrugComb_drugCombCat, DrugBank_ddi, 
                                         by.x = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"), 
                                         by.y = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"), 
                                         all.x = TRUE)


breast_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "BREAST")
lung_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "LUNG")
kidney_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "KIDNEY")
prostate_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "PROSTATE")
ovarian_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "OVARY")
skin_tissue <- which(FIMM_DrugBank_drugCombs_all$tissue == "SKIN")


# Check the associations with different drug combination categories in different cancer types

print(paste0("Number of drug combinations with complete data: ", 
             length(which(complete.cases(FIMM_DrugBank_drugCombs_all)))))

print("Association between total synergistic score and therapeutic efficacy based categories:")
table(FIMM_DrugBank_drugCombs_all$totSyn, FIMM_DrugBank_drugCombs_all$effVsNeff)


print("Association between total synergistic score and pharmacokinetic based categories:")
table(FIMM_DrugBank_drugCombs_all$totSyn, FIMM_DrugBank_drugCombs_all$phamk)

print("Association between total synergistic score and ADV info:")
print("<breast>")
table(FIMM_DrugBank_drugCombs_all$totSyn[breast_tissue], FIMM_DrugBank_drugCombs_all$phamk[breast_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[breast_tissue], FIMM_DrugBank_drugCombs_all$ADR_BreastCancer[breast_tissue])

print("<lung>")
table(FIMM_DrugBank_drugCombs_all$totSyn[lung_tissue], FIMM_DrugBank_drugCombs_all$phamk[lung_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[lung_tissue], FIMM_DrugBank_drugCombs_all$ADR_LungCancer[lung_tissue])

print("<kidney>")
table(FIMM_DrugBank_drugCombs_all$totSyn[kidney_tissue], FIMM_DrugBank_drugCombs_all$phamk[kidney_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[kidney_tissue], FIMM_DrugBank_drugCombs_all$ADR_KidneyCancer[kidney_tissue])

print("<prostate>")
table(FIMM_DrugBank_drugCombs_all$totSyn[prostate_tissue], FIMM_DrugBank_drugCombs_all$phamk[prostate_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[prostate_tissue], FIMM_DrugBank_drugCombs_all$ADR_ProstateCancer[prostate_tissue])
print("<ovarian>")
table(FIMM_DrugBank_drugCombs_all$totSyn[ovarian_tissue], FIMM_DrugBank_drugCombs_all$phamk[ovarian_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[ovarian_tissue], FIMM_DrugBank_drugCombs_all$ADR_OvaryCancer[ovarian_tissue])

print("<skin>")
table(FIMM_DrugBank_drugCombs_all$totSyn[skin_tissue], FIMM_DrugBank_drugCombs_all$phamk[skin_tissue])
table(FIMM_DrugBank_drugCombs_all$totSyn[skin_tissue], FIMM_DrugBank_drugCombs_all$ADR_SkinCancer[skin_tissue])



# Create the cancer type specific drug combinations
# BreastCancer specific
breast_dp = FIMM_DrugBank_drugCombs_all[breast_tissue,]

# LungCancer specific
lung_dp = FIMM_DrugBank_drugCombs_all[lung_tissue,]

# ProstateCancer specific
prostate_dp = FIMM_DrugBank_drugCombs_all[prostate_tissue,]

# KidneyCancer specific
kidney_dp = FIMM_DrugBank_drugCombs_all[kidney_tissue,]

# OvarianCancer specific
ovarian_dp = FIMM_DrugBank_drugCombs_all[ovarian_tissue,]

# SkinCancer specific
skin_dp = FIMM_DrugBank_drugCombs_all[skin_tissue,]


# Save the identified drug pairs
drug_pairs <- list("BreastCancer" = breast_dp,
                   "KidneyCancer" = kidney_dp,
                   "LungCancer" = lung_dp,
                   "OvarianCancer" = ovarian_dp,
                   "ProstateCancer" = prostate_dp,
                   "SkinCancer" = skin_dp)


if(!dir.exists("InputFiles/Drug_combinations/")){
  dir.create("InputFiles/Drug_combinations/", recursive = TRUE)
} 
saveRDS(drug_pairs, "InputFiles/Drug_combinations/Drug_combinations.rds")



print(warnings())