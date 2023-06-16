set.seed(5081)



# dNET RWR on Drug Combinations



# Load libraries
library(optparse)
library(foreach)
library(doParallel)
source("Scripts/Functions/Functions_dNet_RWR_analysis.R")
source("Scripts/Functions/Functions_parallelprocesses.R")





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--nproc"), type = "numeric", default = NULL, 
          help = "Number of processes to use. Default: NULL", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!is.null(opt$nproc)){
  if(!is.numeric(opt$nproc) | (opt$nproc %% 1 != 0)){
    print_help(opt_parser)
    stop("--nproc should be be an integer.", call.=FALSE)
  }
}

# Define global options for this script 
disease <- opt$disease
nproc <- opt$nproc

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Read the drug combinations
drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", disease, ".rds"))
cat("\n\nNumber of drug combinations:\n")
lapply(drugCombs, nrow)



# Read drug target interactions
drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")



# Read the network on which to run RWR
rwr_input_network <- readRDS("InputFiles/Networks/STRING_PPI_Net_database.rds")




# Execute RWR
if(is.null(nproc)){nproc <- detectCores()/2} #Check the number of cores. Use only half of available
cl <- makeCluster(nproc)
registerDoParallel(cl) 


drugCombs_rwr_res1 <- foreach (i=1:nrow(drugCombs$effectiveCombinations),
                               .packages = c("igraph", "tidyr", "dnet")) %dopar% {
                                 tmp <- func_dNetRWR_on_drugCombination(rwr_input_network = rwr_input_network,
                                                                        drug1 = drugCombs$effectiveCombinations[i,c("Drug1_DrugBank_drug_id")], 
                                                                        drug2 = drugCombs$effectiveCombinations[i,c("Drug2_DrugBank_drug_id")],
                                                                        drug_target_ixn = drug_target_ixn,
                                                                        rwr_norm_input = "laplacian",
                                                                        rwr_norm_output = "none", 
                                                                        rwr_restart = 0.5)
                                 tmp
                               }

for(i in 1:length(drugCombs_rwr_res1)){
  names(drugCombs_rwr_res1)[i] <- paste("Eff", drugCombs$effectiveCombinations[i,c("Drug1_DrugBank_drug_id")], 
                                        drugCombs$effectiveCombinations[i,c("Drug2_DrugBank_drug_id")], sep = "__")
}



drugCombs_rwr_res2 <- foreach (i=1:nrow(drugCombs$adverseCombinations),
                               .packages = c("igraph", "tidyr", "dnet")) %dopar% {
                                 tmp <- func_dNetRWR_on_drugCombination(rwr_input_network = rwr_input_network,
                                                                        drug1 = drugCombs$adverseCombinations[i,c("Drug1_DrugBank_drug_id")], 
                                                                        drug2 = drugCombs$adverseCombinations[i,c("Drug2_DrugBank_drug_id")],
                                                                        drug_target_ixn = drug_target_ixn,
                                                                        rwr_norm_input = "laplacian",
                                                                        rwr_norm_output = "none", 
                                                                        rwr_restart = 0.5)
                                 tmp
                               }

for(i in 1:length(drugCombs_rwr_res2)){
  names(drugCombs_rwr_res2)[i] <- paste("Adv", drugCombs$adverseCombinations[i,c("Drug1_DrugBank_drug_id")], 
                                        drugCombs$adverseCombinations[i,c("Drug2_DrugBank_drug_id")], sep = "__")
}


stopCluster(cl)
unregister_dopar()





# Export results
drugCombs_rwr_res_final <- list()
drugCombs_rwr_res_final$effectiveCombinations <- drugCombs_rwr_res1
drugCombs_rwr_res_final$adverseCombinations <- drugCombs_rwr_res2

if(!dir.exists(paste0("OutputFiles/Model_train/", disease, "/"))){
  dir.create(paste0("OutputFiles/Model_train/", disease, "/"), recursive = TRUE)
} 

save(drugCombs_rwr_res_final, file = paste0("OutputFiles/Model_train/", disease, "/dNetRWR050_", disease, ".rda"))



print(warnings())