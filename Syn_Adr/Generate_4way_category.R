set.seed(5081)



# Script to generate effective vs adverse categorisation for drug combinations



# Load libraries
library(unixtools)
library(optparse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


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
cat(paste0("\nDisease: ", disease))


# Read the synergy level of the  drug combination
drugCombs_data <- readRDS(paste0("InputFiles/Drug_combination_data/drugCombs_data_", disease, ".rds"))
drugCombs_data <- drugCombs_data[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "Syn_level")]


# Read the DDI data
DrugBank_ddi <- readRDS("InputFiles/Reference_list/DrugBank_DDI_processed.rds")
DrugBank_ddi <- DrugBank_ddi[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", paste0("ADR_", disease))]


# Merge the data
drugCombs_cat <- merge(DrugBank_ddi, drugCombs_data, 
                       by = c("Drug1_DrugBank_id", "Drug2_DrugBank_id"),
                       all.y = TRUE)


# Assign the drug combination categories
drugCombs_cat$class_4way <- NA

drugCombs_cat[drugCombs_cat$Syn_level >= 3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown",]$class_4way <- "Syn_Eff"
drugCombs_cat[drugCombs_cat$Syn_level >= 3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive",]$class_4way <- "Syn_Adv"

drugCombs_cat[drugCombs_cat$Syn_level <= -3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "unknown",]$class_4way  <- "Ant_Eff"
drugCombs_cat[drugCombs_cat$Syn_level <= -3 & drugCombs_cat[, paste0("ADR_", disease)] %in% "adr_positive",]$class_4way  <- "Ant_Adv"




# drugCombs_cat <- drugCombs_cat[!is.na(drugCombs_cat$class_4way), c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_4way")]
drugCombs_cat <- drugCombs_cat[, c("Drug1_DrugBank_id", "Drug2_DrugBank_id", "class_4way")]


if(!dir.exists("Syn_Adr/Drug_combination_class/")){
  dir.create("Syn_Adr/Drug_combination_class/", recursive = TRUE)
}
saveRDS(drugCombs_cat, file = paste0("Syn_Adr/Drug_combination_class/drugCombs_cat_4way_", disease, ".rds"))


print(warnings())