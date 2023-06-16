set.seed(5081)



# Perform FGSEA on drug combinations based on the RWR results



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("Scripts/Functions/Functions_parallelprocesses.R")





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
          help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer.", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Load RWR results
load(file = paste0("OutputFiles/Model_train/", disease, "/dNetRWR050_", disease, ".rda"))


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

NES_WithdrawalAdr2Gene <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_WithdrawalAdr2Gene)[1] <- "ADR"





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

NES_CombinedDisAdr2Gene <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_CombinedDisAdr2Gene)[1] <- "Term"



# Final results
fgsea_result <- list()
fgsea_result$NES_Disease2Gene <- NES_Disease2Gene
fgsea_result$NES_WithdrawalAdr2Gene <- NES_WithdrawalAdr2Gene
fgsea_result$NES_CombinedDisAdr2Gene <- NES_CombinedDisAdr2Gene

saveRDS(fgsea_result, file = paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_EfficacySafety_", disease, ".rds"))



print(warning())