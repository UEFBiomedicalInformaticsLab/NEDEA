set.seed(5081)
rm(list = ls())



# Script to split data into training and test set for ML (version 5)


# Load libraries
library(optparse)
library(caret)




# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character"),
  make_option(c("--train_perc"), type = "numeric", default = 70, 
              help = "Percentage of the data to be used as training data. Should be an integer between 1 and 100.", metavar = "numeric"),
  make_option(c("--splits"), type = "numeric", default = 10, 
              help = "Number of training/test data splits to me made. Should be an integer.", metavar = "numeric")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}


if(!is.numeric(opt$train_perc) | (opt$train_perc %% 1 != 0)){
  print_help(opt_parser)
  stop("--train_perc should be be an integer between 1 and 100.", call.=FALSE)
}


if(!is.numeric(opt$splits) | (opt$splits %% 1 != 0)){
  print_help(opt_parser)
  stop("--splits should be be an integer.", call.=FALSE)
}



# Define global options for this script 
disease <- opt$disease
train_perc <- opt$train_perc
splits <- opt$splits

cat(paste0("\n\nExecuting for: ", disease, "\n\n"))



# Read the drug combinations
drugCombs <- readRDS(paste0("InputFiles/DrugCombinations/DrugCombs_v5/DrugComb_", disease, "_v5.rds"))
cat("\n\nNumber of drug combinations:\n")
lapply(drugCombs, nrow)


drugCombs <- do.call(rbind, drugCombs)
drugCombs$Class <- substr(x = row.names(drugCombs), start = 1, stop = 3) 
drugCombs$Class <- gsub(pattern = "^eff", replacement = "Eff", x = drugCombs$Class)
drugCombs$Class <- gsub(pattern = "^adv", replacement = "Adv", x = drugCombs$Class)
drugCombs$Name <- paste(drugCombs$Class, drugCombs$Drug1_DrugBank_drug_id, drugCombs$Drug2_DrugBank_drug_id, sep = "__")
row.names(drugCombs) <- NULL


trainIndex <- createDataPartition(y = drugCombs$Class  , p = train_perc/100, list = FALSE, times = splits)


ML_data_split <- list()
for(i in colnames(trainIndex)){
  trainData <- drugCombs[trainIndex[, i], ]
  testData <- drugCombs[-trainIndex[, i], ]
  
  ML_data_split[[i]] <- list("train" = trainData, "test" = testData)
  }

saveRDS(ML_data_split, file = paste0("Analysis/STRING/DrugCombs_v5/", disease, "/ML_dataSplit_", disease, ".rds"))
