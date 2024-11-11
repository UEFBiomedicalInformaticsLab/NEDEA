set.seed(5081)


# Custom script to execute SAveRUNNER
# Notes:
# (a) Reads the input files, sets config details and executes for each disease


#####


main_exc_dir <- getwd()

for(disease in c("BreastCancer", "KidneyCancer", "LungCancer", "OvaryCancer", "ProstateCancer", "SkinCancer")){
  
  
  print(paste0(disease, "-------------"))
  
  # Prepare the input files 
  
  #interactome
  interactome <- read.table("InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_Interactome.tsv", header = T, sep = '\t', check.names = F, quote = "")
  
  # disease gene
  disease_gene <- read.table(paste0("InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_PhenotypeGenes_", disease, ".tsv"), header = T, sep = '\t', check.names = F, quote = "")
  disease_gene <- unique(disease_gene)
  disease_gene$disease <- gsub(pattern = "\"", x = disease_gene$disease, replacement = "")
  
  # drug target 
  drug_target <- read.table(paste0("InputFiles/Compare_external/SAveRUNNER/SAveRUNNER_DrugCombTargets_", disease, ".tsv"), header = T, sep = '\t', check.names = F, quote = "")
  drug_target <- unique(drug_target)
  drug_target$Drug <- tolower(drug_target$Drug)
  
  input_file <-list(interactome = interactome,
                    disease_gene = disease_gene,
                    drug_target = drug_target)
  
  
  #####
  
  
  # Prepare the config file
  
  
  # parameters for computing start network with edge-weight = proximity
  diseases <- unique(disease_gene$disease)
  
  # parameters for computing end network
  dirRes <- paste0(main_exc_dir, "/OutputFiles/Compare_external/SAveRUNNER/", disease, "/")       
  
  interaction = "similarity"  # edge-weight = similarity or proximity
  pval_thr = 0.05             # select significative drug-disease association
  adjust_link = T            # adjust similarity or not
  new_link = F                # add new drug-disease association or not (without compute pval)
  
  # parameters for making figure
  if( (interaction == "proximity") ) distance = "proximity"
  if( (interaction == "similarity") & (adjust_link == F) ) distance = "similarity"
  if( (interaction == "similarity") & (adjust_link == T) ) distance = "adjusted_similarity"
  
  # parameters for computing subnetwork
  sel_drug = NULL
  sel_disease = NULL
  
  input_parameter <- list(diseases = diseases,
                          dirRes = dirRes,
                          interaction = interaction,
                          pval_thr = pval_thr, 
                          adjust_link = adjust_link,
                          new_link = new_link,
                          distance = distance,
                          sel_drug = sel_drug,
                          sel_disease = sel_disease)
  
  
  #####
  
  
  # Prepare the main 
  
  options(stringsAsFactors = F)
  
  setwd("External_tools/SAveRUNNER/code/")
  
  source("src/script/getLibrary.R")
  source("src/script/getSource.R")
  
  getLibrary()
  getSource()
  
  # -- AG
  # Already defined above
  # In source the two functions have been removed.
  # input_parameter <- config()
  # input_file <- inputFiles()
  # -- AG
  
  
  # -- AG
  if(!dir.exists(input_parameter$dirRes)){
    dir.create(input_parameter$dirRes, recursive = TRUE)
  } 
  # -- AG
  
  output_file <- networkFiles()
  
  # compute drug-disease network
  destfile = output_file$filename_out_allPval
  if ( !file.exists(destfile) ) mainStartNetwork()
  
  # select significative drug-disease association and (or not) adjust them
  mainEndNetwork()
  
  # make figure
  
  if( length(input_parameter$diseases) > 1) mainFigure()
  
  # create disease specific subnetwork
  if( !is.null(input_parameter$sel_disease) ) mainSubnetwork()
  
  
  setwd(main_exc_dir)
  
  file.copy("External_tools/SAveRUNNER/code/Drug_Disease_network.txt", paste0("OutputFiles/Compare_external/SAveRUNNER/", disease, "/"))  
  file.remove("External_tools/SAveRUNNER/code/Drug_Disease_network.txt")
  
}


#####


print(warnings())