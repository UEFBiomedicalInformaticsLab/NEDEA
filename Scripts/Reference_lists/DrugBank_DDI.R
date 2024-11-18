set.seed(5081)


# Script to analyse the drug-drug interactions (DDIs) reported in DrugBank


# Load libraries
library(tidyverse)


#####


# Read the DrugBank data
DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")


# Read the drug type information
DrugBank_drug_type <- DrugBank_data$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


# Read the DDI data
DrugBank_ddi <- DrugBank_data$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% DrugBank_drug_type$primary_key & DrugBank_ddi$Drug2_DrugBank_id %in% DrugBank_drug_type$primary_key, ] # retain only small molecular drugs


#####


# Mark which drug combinations have serious ADR
DrugBank_ddi$ADR_status <- "adr_negative"
DrugBank_ddi[grepl("risk or severity of .+toxicity can be increased|risk or severity of liver damage can be increased|risk or severity of adverse effects can be increased |increase the .+toxic activities", DrugBank_ddi$description, ignore.case = TRUE), "ADR_status"] <- "adr_positive"


#####


if(!dir.exists("InputFiles/Reference_list/"))dir.create("InputFiles/Reference_list/", recursive = TRUE)
saveRDS(DrugBank_ddi, "InputFiles/Reference_list/DrugBank_DDI_processed.rds")


#####


print(warnings())