set.seed(5081)


start_time <- Sys.time()
print(paste0("Start time: ", start_time))

# Load libraries
library(optparse)
source("Scripts/Functions/Functions_FGSEA.R")



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--drug_target_type"), type = "character", default = "known", 
              help = "The type of drug target to use. Possible values: known, PS, SIGNOR, NPA, RI, KEGG, all. Default: known", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call.=FALSE)
}


# Define global options for this script 
disease <- opt$disease
drug_target_type <- opt$drug_target_type


cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease))
cat(paste0("\nDrug target type: ", drug_target_type, "\n\n"))


# Read the RWR results
rwr_result_file_path <- paste0("OutputFiles/RWR_results/rwrProbs_", disease, "_", drug_target_type, ".rds")
if(!file.exists(rwr_result_file_path)){
  stop(paste0("Missing file. Check if \'", drug_target_type, "\' was used to compile RWR."), call. = TRUE)
}
rwr_result <- readRDS(file = rwr_result_file_path)


fgsea_result_final <- list()

# FGSEA on efficacy library
cat("\n--- Executing FGSEA on efficacy library")
enrichment_lib <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib) <- paste0("[DISEASE] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib)

fgsea_result_final[["efficacy"]] <- fgsea_result


# FGSEA on safety library
cat("\n--- Executing FGSEA on safety library")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/drugWithdrawal_Adr2Gene_lib.rds")
names(enrichment_lib) <- paste0("[ADR] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib)

fgsea_result_final[["safety"]] <- fgsea_result


# FGSEA on combined efficacy-safety library
cat("\n--- Executing FGSEA on combined efficacy-safety library")
enrichment_lib_1 <- readRDS(paste0("InputFiles/Enrichment_analysis_libraries/Disease2Gene_", disease, "_lib.rds"))
names(enrichment_lib_1) <- paste0("[DISEASE] ", names(enrichment_lib_1))
enrichment_lib_2 <- readRDS("InputFiles/Enrichment_analysis_libraries/drugWithdrawal_Adr2Gene_lib.rds")
names(enrichment_lib_2) <- paste0("[ADR] ", names(enrichment_lib_2))
enrichment_lib <- c(enrichment_lib_1, enrichment_lib_2)

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib)

fgsea_result_final[["combinedEfficacySafety"]] <- fgsea_result




if(!dir.exists("OutputFiles/FGSEA_results/")){
  dir.create("OutputFiles/FGSEA_results/", recursive = TRUE)
} 
saveRDS(fgsea_result_final, paste0("OutputFiles/FGSEA_results/fgseaNES_EfficacySafety_", disease, "_", drug_target_type, ".rds"))


end_time <- Sys.time()
print(paste0("Start time: ", end_time))
print(paste0("Time taken: ", round(end_time - start_time, 2)))

print(warnings())