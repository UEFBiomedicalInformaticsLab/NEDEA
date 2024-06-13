
disease <- "SkinCancer"







# Read the DDI data
DrugBank_ddi <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_ddi <- DrugBank_ddi$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")


# Read the drug combination category
drugCombs_cat <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs_cat$comb_name <- paste(drugCombs_cat$Drug1_DrugBank_id, drugCombs_cat$Drug2_DrugBank_id, sep = "_")
# drugCombs_cat <- drugCombs_cat[, c("comb_name", "class_EffAdv")]
drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_EffAdv), ]
row.names(drugCombs_cat) <- NULL



drugCombs_cat$DDI <- c()
for(i in 1:nrow(drugCombs_cat)){
  drug1 <- drugCombs_cat[i, "Drug1_DrugBank_id"]
  drug2 <- drugCombs_cat[i, "Drug2_DrugBank_id"]
  
  tmp1 <- unique( DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% c(drug1, drug2) & DrugBank_ddi$Drug2_DrugBank_id %in% c(drug1, drug2), "description"] )
  drugCombs_cat[i, "DDI"] <- paste(tmp1, collapse = "; ")
}


View(drugCombs_cat)

