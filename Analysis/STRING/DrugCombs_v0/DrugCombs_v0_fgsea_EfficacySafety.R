set.seed(5081)
rm(list = ls())



# Perform FGSEA on drug combinations based on the RWR results (version 0)



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("/research/groups/fortino/arindam/DrugCombination_1/Scripts/Functions/Functions_parallelprocesses.R")





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

# Define global options for this script 
disease <- opt$disease
disease <- "LungCancer"
cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Load RWR results
load(file = paste0("Analysis/STRING/DrugCombs_v0/", disease, "/dNetRWR050_", disease, ".rda"))


effectiveComb_rwr <- drugCombs_rwr_res_final$effectiveCombinations
adverseComb_rwr <- drugCombs_rwr_res_final$adverseCombinations

# Check if all resuls have valied number of columns
effectiveComb_rwr <- effectiveComb_rwr[lapply(effectiveComb_rwr, ncol) == 7]
adverseComb_rwr <- adverseComb_rwr[lapply(adverseComb_rwr, ncol) == 7]
rm(drugCombs_rwr_res_final)


cat("\n\nNumbe rof drug combinations: ")
cat(paste0("\n - Effective combinations = ", length(effectiveComb_rwr)))
cat(paste0("\n - Adverse combinations = ", length(adverseComb_rwr), "\n\n"))





# Run FGSEA on Disease2Gene library (efficacy)
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))


effectiveComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                        rwr_data = effectiveComb_rwr, 
                                                        quantile_prob = 0.9)

adverseComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                      rwr_data = adverseComb_rwr, 
                                                      quantile_prob = 0.9)


effectiveComb_NES <- func_extract_fgsea_result(enrichment_result = effectiveComb_enrichment,
                                               result_type = "NES",
                                               enrichment_library = enrichment_lib)

adverseComb_NES <- func_extract_fgsea_result(enrichment_result = adverseComb_enrichment,
                                             result_type = "NES",
                                             enrichment_library = enrichment_lib)

NES_Disease2Gene <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_Disease2Gene)[1] <- "Disease"





# Run FGSEA on drug withdrawal based Adr2Gene library (safety)
enrichment_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")


effectiveComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                        rwr_data = effectiveComb_rwr, 
                                                        quantile_prob = 0.9)

adverseComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                      rwr_data = adverseComb_rwr, 
                                                      quantile_prob = 0.9)


effectiveComb_NES <- func_extract_fgsea_result(enrichment_result = effectiveComb_enrichment,
                                               result_type = "NES",
                                               enrichment_library = enrichment_lib)

adverseComb_NES <- func_extract_fgsea_result(enrichment_result = adverseComb_enrichment,
                                             result_type = "NES",
                                             enrichment_library = enrichment_lib)

NES_withdrawalAdr2Gene <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_withdrawalAdr2Gene)[1] <- "ADR"





# Run FGSEA on combined Disease2Gene and drug withdrawal based Adr2Gene library (efficacy + safety)
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))

enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)


effectiveComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                        rwr_data = effectiveComb_rwr, 
                                                        quantile_prob = 0.9)

adverseComb_enrichment <- func_fgsea_from_rwr_probCut(enrichment_library = enrichment_lib, 
                                                      rwr_data = adverseComb_rwr, 
                                                      quantile_prob = 0.9)


effectiveComb_NES <- func_extract_fgsea_result(enrichment_result = effectiveComb_enrichment,
                                               result_type = "NES",
                                               enrichment_library = enrichment_lib)

adverseComb_NES <- func_extract_fgsea_result(enrichment_result = adverseComb_enrichment,
                                             result_type = "NES",
                                             enrichment_library = enrichment_lib)

DiseaseAdr2Gene <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(DiseaseAdr2Gene)[1] <- "Term"





# Final results
fgsea_result <- list()
fgsea_result$NES_Disease2Gene <- NES_Disease2Gene
fgsea_result$NES_withdrawalAdr2Gene <- NES_withdrawalAdr2Gene
fgsea_result$DiseaseAdr2Gene <- DiseaseAdr2Gene

saveRDS(fgsea_result, file = paste0("Analysis/STRING/DrugCombs_v0/", disease, "/fgseaProbCut_", disease, ".rds"))

print(warning())