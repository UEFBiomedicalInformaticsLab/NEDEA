set.seed(5081)


# Perform FGSEA on drug combinations based on the RWR results (version 5)





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





# Run FGSEA on CHG_keggPath2Gene_lib
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/CHG_keggPath2Gene_lib.rds"))
cat("\n\n\n- Executing FGSEA for CHG_keggPath2Gene_lib\n")


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

NES_keggPath <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_keggPath)[1] <- "Pathway"





# Run FGSEA on msigdb_ReactomePath2Gene_lib
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/msigdb_ReactomePath2Gene_lib.rds"))
cat("\n\n\n- Executing FGSEA for msigdb_ReactomePath2Gene_lib\n")


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

NES_ReactomePath <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_ReactomePath)[1] <- "Pathway"






# Run FGSEA on SMPDb_Pathway2Gene_lib (Drug Metabolism)
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Metabolism`
cat("\n\n\n- Executing FGSEA for SMPDb_Pathway2Gene_lib (Drug metabolism)\n")


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

NES_SMPDbPath_DrugMet <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_SMPDbPath_DrugMet)[1] <- "Pathway"





# Run FGSEA on SMPDb_Pathway2Gene_lib (Drug Action)
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/SMPDb_Pathway2Gene_lib.rds"))
enrichment_lib <- enrichment_lib$`Drug Action`
cat("\n\n\n- Executing FGSEA for SMPDb_Pathway2Gene_lib (Drug Action)\n")


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

NES_SMPDbPath_DrugAction <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_SMPDbPath_DrugAction)[1] <- "Pathway"





# Run FGSEA on miscellaneous_gene_lib
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_Analysis_Libraries/miscellaneous_gene_lib.rds"))
cat("\n\n\n- Executing FGSEA for miscellaneous_gene_lib\n")


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

NES_miscGeneSet <- merge(effectiveComb_NES, adverseComb_NES, by = 0 , all = TRUE)
names(NES_miscGeneSet)[1] <- "GeneSet"





# Final results
fgsea_result <- list()
fgsea_result$NES_keggPath <- NES_keggPath
fgsea_result$NES_ReactomePath <- NES_ReactomePath
fgsea_result$NES_SMPDbPath_DrugMet <- NES_SMPDbPath_DrugMet
fgsea_result$NES_SMPDbPath_DrugAction <- NES_SMPDbPath_DrugAction
fgsea_result$NES_miscGeneSet <- NES_miscGeneSet

saveRDS(fgsea_result, file = paste0("OutputFiles/Model_train/", disease, "/fgseaProbCut_CommonLib_", disease, ".rds"))



print(warning())