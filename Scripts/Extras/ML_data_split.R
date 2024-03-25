set.seed(5081)



# Script to split data into training and test set for ML



# Load libraries
library(optparse)
library(caret)



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--folds"), type = "numeric", default = 5, 
              help = "Number of folds to be created. Should be an integer.", metavar = "numeric"),
  make_option(c("--repeats"), type = "numeric", default = 5, 
              help = "Number of repeats of fold splitting to be performed. Should be an integer.", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

if(!is.numeric(opt$folds) | (opt$folds %% 1 != 0)){
  print_help(opt_parser)
  stop("--folds should be be an integer.", call.=FALSE)
}

if(!is.numeric(opt$repeats) | (opt$repeats %% 1 != 0)){
  print_help(opt_parser)
  stop("--repeats should be be an integer.", call.=FALSE)
}



# Define global options for this script 
disease <- opt$disease
folds <- opt$folds
repeats <- opt$repeats

cat(paste0("\n\nExecuting for: ", disease, "\n"))
cat(paste0("Folds: ", folds, "\n"))
cat(paste0("Repeats: ", repeats, "\n\n"))



# Read the drug combinations
drugCombs <- readRDS(paste0("InputFiles/Drug_combination_class/drugCombs_cat_effVadv_", disease, ".rds"))
drugCombs <- drugCombs[!is.na(drugCombs$class_EffAdv), ]
drugCombs$comb_name <- paste(drugCombs$Drug1_DrugBank_id, drugCombs$Drug2_DrugBank_id, sep = "_")
cat("\n\nNumber of drug combinations:\n")
print(table(drugCombs$class_EffAdv))


# Creating folds with repeats

ML_data_split_repeat <- list()
ML_data_split <- list()

for(j in seq(repeats)){
  
  set.seed(5081 + j)
  
  fldIndex <- createFolds(y = drugCombs$class_EffAdv, k = folds, list = TRUE, returnTrain = FALSE)
  
  for(i in 1:length(fldIndex)) {
    testData <- drugCombs[fldIndex[[i]], ]
    remainingFolds <- fldIndex[-i]
    trainData <- drugCombs[Reduce(union, lapply(remainingFolds, unlist)), ]
    ML_data_split[[i]] <- list("train" = trainData, "test" = testData)
  }
  names(ML_data_split) <- names(fldIndex)
  ML_data_split_repeat[[j]] <- ML_data_split
}
names(ML_data_split_repeat) <- paste0("Repeat", seq(repeats))
ML_data_split_repeat <- unlist(ML_data_split_repeat, recursive = FALSE, )

names(ML_data_split_repeat) <- gsub("\\.", "_", names(ML_data_split_repeat))




if(!dir.exists("OutputFiles/ML_data_split/")){dir.create("OutputFiles/ML_data_split/", recursive = TRUE)}
saveRDS(ML_data_split_repeat, file = paste0("OutputFiles/ML_data_split/ML_dataSplit_", disease, ".rds"))



print(warnings())