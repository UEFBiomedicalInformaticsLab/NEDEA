set.seed(5081)



# Script to execute FGSEA on the RWR results using miscellaneous gene set library


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

# FGSEA on miscellaneous gene set library
cat("\n--- Executing FGSEA on miscellaneous gene set library")
enrichment_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/miscellaneous_gene_lib.rds")
names(enrichment_lib) <- paste0("[MISC.] ", names(enrichment_lib))

fgsea_result <- func_run_FGSEA_on_RWR(rwr_data = rwr_result, 
                                      enrichment_library = enrichment_lib)

fgsea_result_final[["misc"]] <- fgsea_result


if(!dir.exists("OutputFiles/FGSEA_results/")){
  dir.create("OutputFiles/FGSEA_results/", recursive = TRUE)
} 
saveRDS(fgsea_result_final, paste0("OutputFiles/FGSEA_results/fgseaNES_Miscellaneous_", disease, "_", drug_target_type, ".rds"))



print(warnings())