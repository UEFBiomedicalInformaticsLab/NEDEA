set.seed(5081)



# Script to compile results from external predictions using the CDCDB drug combinations (AACT)
# Notes:
# (a) In the final output, checks if the drugs are already approved for certain type of cancers.


# Load libraries
library(openxlsx)
library(tidyverse)
library(sparklyr)
library(sparklyr.nested)



# Compile the results into one file
files <- list.files(path = "OutputFiles/External_predictions/CDCDB_AACT",
                    pattern = "^ExtPred_CDCDB_AACT__", 
                    full.names = TRUE)


extPred_drugCombs <- data.frame()
for(file in files){
  disease <- strsplit(x = file, split = "__")[[1]][[2]]
  disease <- strsplit(x = disease, split = "_")[[1]][[1]]
  
  model <- strsplit(x = file, split = "__")[[1]][[2]]
  model <- strsplit(x = model, split = "_")[[1]][[2]]
  
  balance <- strsplit(x = file, split = "__")[[1]][[2]]
  balance <- strsplit(x = balance, split = "_")[[1]][[3]]
  
  featureType <- strsplit(x = file, split = "__")[[1]][[2]]
  featureType <- strsplit(x = featureType, split = "_")[[1]][[4]]
  featureType <- strsplit(x = featureType, split = "\\.")[[1]][[1]]
  
  
  tmp1 <- read.xlsx(file)
  tmp1$disease <- disease
  tmp1$model <- model
  tmp1$balance <- balance
  tmp1$featureType <- featureType
  
  extPred_drugCombs <- rbind(extPred_drugCombs, tmp1)
}
rm(tmp1)



# Save as file
if(!dir.exists("OutputFiles/External_predictions/CDCDB_AACT/")){
  dir.create("OutputFiles/External_predictions/CDCDB_AACT/", recursive = TRUE)
}                          
write.xlsx(extPred_drugCombs, "OutputFiles/External_predictions/CDCDB_AACT/Compiled_ExtPred_CDCDB_AACT.xlsx", 
           rowNames = FALSE, overwrite = TRUE)



print(warnings())