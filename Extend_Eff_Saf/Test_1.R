
drug_target_type <- "known"
disease <- "BreastCancer"

rwr_data = rwr_result
enrichment_library = enrichment_lib




# Extract the thresholds for RWR probabilities
rwr_threshold <- sapply(apply(rwr_data, 2, func_RWR_threshold), function(x) x$ELB)


drugComb <- "DB04868_DB06813"


rwr_data_select <- rwr_data[, drugComb]
rankedGeneList <- sort(rwr_data_select[rwr_data_select > rwr_threshold[drugComb]], decreasing = TRUE)
rankedGeneList



res <- fgseaMultilevel(pathways = enrichment_library,
                       stats = rankedGeneList,
                       minSize = 5, 
                       maxSize = Inf, 
                       scoreType = "pos", 
                       BPPARAM = SerialParam())

table(names(rankedGeneList) %in% enrichment_library$Efficacy)

table(names(rankedGeneList) %in% enrichment_library$Safety)





fgsea_result_final <- readRDS(paste0("Extend_Eff_Saf/FGSEA_results/fgseaNES_extendedEfficacySafety_", disease, "_", drug_target_type, ".rds"))
fgsea_result_final <- fgsea_result_final$extendedEfficacySafety



tmp1 <- colSums(fgsea_result_final)
tmp1[tmp1 > 0]

fgsea_result_final <- fgsea_result_final[, drugComb, drop = FALSE]


