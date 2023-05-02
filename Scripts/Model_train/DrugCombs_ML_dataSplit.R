# Script to split data into training and test set for ML



# Load libraries
library(optparse)
library(caret)





# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--folds"), type = "numeric", default = 5, 
              help = "Number of folds to be created. Should be an integer.", metavar = "numeric")
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



# Define global options for this script 
disease <- opt$disease
folds <- opt$folds

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Read the drug combinations
drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugComb_", disease, ".rds"))
cat("\n\nNumber of drug combinations:\n")
lapply(drugCombs, nrow)


drugCombs <- do.call(rbind, drugCombs)
drugCombs$Class <- substr(x = row.names(drugCombs), start = 1, stop = 3) 
drugCombs$Class <- gsub(pattern = "^eff", replacement = "Eff", x = drugCombs$Class)
drugCombs$Class <- gsub(pattern = "^adv", replacement = "Adv", x = drugCombs$Class)
drugCombs$Name <- paste(drugCombs$Class, drugCombs$Drug1_DrugBank_drug_id, drugCombs$Drug2_DrugBank_drug_id, sep = "__")
row.names(drugCombs) <- NULL



# Creating folds
fldIndex <- createFolds(y = drugCombs$Class, k = folds, list = TRUE, returnTrain = FALSE)

ML_data_split <- list()
for(i in 1:length(fldIndex)) {
  testData <- drugCombs[fldIndex[[i]], ]
  remainingFolds <- fldIndex[-i]
  trainData <- drugCombs[Reduce(union, lapply(remainingFolds, unlist)), ]
  ML_data_split[[i]] <- list("train" = trainData, "test" = testData)
}
names(ML_data_split) <- names(fldIndex)

saveRDS(ML_data_split, file = paste0("OutputFiles/Model_train/", disease, "/ML_dataSplit_", disease, ".rds"))



print(warnings())