set.seed(5081)


# Custom script to execute SAveRUNNER

# Load libraries
# library(reticulate)

#####


main_exc_dir <- getwd()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  print(paste0(disease, "-------------"))
  
  if(!dir.exists(paste0("OutputFiles/Compare_external/deepDDI/", disease))){
    dir.create(paste0("OutputFiles/Compare_external/deepDDI/", disease), recursive = TRUE)
  } 
  
  setwd("External_tools/deepddi/")
  
  
  cmd <- paste0("python run_DeepDDI.py -i ", main_exc_dir, "/InputFiles/Compare_external/deepDDI/deepDDI_DrugCombSMILES_", disease, ".txt -o ", main_exc_dir, "/OutputFiles/Compare_external/deepDDI/", disease)
  
  system(cmd)
  
  setwd(main_exc_dir)
}


#####


print(warnings())