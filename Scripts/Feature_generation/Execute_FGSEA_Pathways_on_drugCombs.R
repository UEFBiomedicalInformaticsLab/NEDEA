set.seed(5081)



# Script to execute FGSEA on the RWR results using pathway library


# Load libraries
library(unixtools)
library(optparse)
source("Scripts/Functions/Functions_FGSEA.R")

# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


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
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$drug_target_type %in% c("known", "PS", "SIGNOR", "NPA", "RI", "KEGG", "all")){
  print_help(opt_parser)
  stop("--drug_target_type should be: known, PS, SIGNOR, NPA, RI, KEGG, all", call. = FALSE)
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

# FGSEA on KEGG pathways liked to cancer hallmarks
cat("\n--- Executing FGSEA on KEGG pathways liked to cancer hallmarks")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/CHG_keggPath2Gene_lib.rds")
names(enrichment_lib) <- paste0("[KEGG] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib,
                                      disease = disease,
                                      drug_target_type = drug_target_type)

fgsea_result_final[["kegg"]] <- fgsea_result


# FGSEA on SMPDB pathways (drug metabolism)
cat("\n--- Executing FGSEA on SMPDB pathways (drug metabolism)")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds")
enrichment_lib <- enrichment_lib$`Drug Metabolism`
names(enrichment_lib) <- paste0("[SMPDB_DRUG_METABOLISM] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib,
                                      disease = disease,
                                      drug_target_type = drug_target_type)

fgsea_result_final[["smpdbDrugMet"]] <- fgsea_result


# FGSEA on SMPDB pathways (drug action)
cat("\n--- Executing FGSEA on SMPDB pathways (drug action)")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/SMPDb_Pathway2Gene_lib.rds")
enrichment_lib <- enrichment_lib$`Drug Action`
names(enrichment_lib) <- paste0("[SMPDB_DRUG_ACTION] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib,
                                      disease = disease,
                                      drug_target_type = drug_target_type)

fgsea_result_final[["smpdbDrugAct"]] <- fgsea_result




if(!dir.exists("OutputFiles/FGSEA_results/")){
  dir.create("OutputFiles/FGSEA_results/", recursive = TRUE)
} 
saveRDS(fgsea_result_final, paste0("OutputFiles/FGSEA_results/fgseaNES_Pathway_", disease, "_", drug_target_type, ".rds"))



print(warnings())