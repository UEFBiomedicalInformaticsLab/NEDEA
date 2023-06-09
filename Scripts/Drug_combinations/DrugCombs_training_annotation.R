# Script to create annotations of training drug combinations

library(tidyverse)
library(openxlsx)

drugCombs <- list()
drugCombs[["LungCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_LungCancer.rds")
drugCombs[["BreastCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_BreastCancer.rds")
drugCombs[["KidneyCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_KidneyCancer.rds")
drugCombs[["OvaryCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_OvaryCancer.rds")
drugCombs[["ProstateCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_ProstateCancer.rds")
drugCombs[["SkinCancer"]] <- readRDS("InputFiles/DrugCombinations/DrugComb_SkinCancer.rds")

drugCombs <- lapply(drugCombs, function(x){bind_rows(x = x, .id = "Class")})

DrugBank_Drugs <- read.csv("Databases/DrugBank/drug.csv", header = TRUE)

drugCombs <- lapply(drugCombs, function(x){
    x$Drug1_name <- DrugBank_Drugs$name[match(x$Drug1_DrugBank_drug_id, DrugBank_Drugs$primary_key)]
    x$Drug2_name <- DrugBank_Drugs$name[match(x$Drug2_DrugBank_drug_id, DrugBank_Drugs$primary_key)]
    x$Class <- gsub("effectiveCombinations", "Eff", x$Class)
    x$Class <- gsub("adverseCombinations", "Adv", x$Class)
    x
})
write.xlsx(drugCombs, "InputFiles/DrugCombinations/DrugCombs_training_annotation.xlsx", overwrite = TRUE)