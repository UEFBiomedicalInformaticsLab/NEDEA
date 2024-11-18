set.seed(5081)


# Script to generate effective vs adverse categorisation for drug combinations


# Load libraries
library(unixtools)
library(optparse)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL,
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)


opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call. = FALSE)
}

if(!opt$disease %in% c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  print_help(opt_parser)
  stop("--disease must be one of the following: BreastCancer, KidneyCancer, LungCancer, OvaryCancer, ProstateCancer, SkinCancer", call. = FALSE)
}


# Define global options for this script
disease <- opt$disease

cat("\n\nUsing the following parameters: ")
cat(paste0("\nDisease: ", disease, "\n"))


#####


# Read the synergy level of the  drug combination
drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level")]


# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
DrugBank_ddi <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "ADR_status")]


# Merge the data
drugCombs_cat <- drugCombs_data %>% 
  left_join(DrugBank_ddi, c("Drug1_DrugBank_id", "Drug2_DrugBank_id")) %>% 
  as.data.frame()


# Assign the drug combination categories
drugCombs_cat$class_EffAdv <- NA
drugCombs_cat[drugCombs_cat$Syn_level >= 3 & drugCombs_cat$ADR_status %in% "adr_negative",]$class_EffAdv <- "Eff"
drugCombs_cat[drugCombs_cat$Syn_level <= -2 & drugCombs_cat$ADR_status %in% "adr_positive",]$class_EffAdv  <- "Adv"
drugCombs_cat <- drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_EffAdv")]


if(!dir.exists("InputFiles/Drug_combination_class/")){
  dir.create("InputFiles/Drug_combination_class/", recursive = TRUE)
}
saveRDS(drugCombs_cat, file = paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))


#####


print(warnings())