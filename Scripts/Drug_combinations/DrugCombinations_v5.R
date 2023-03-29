set.seed(5081)
rm(list = ls())



# Drug combinations (version 5)



# Notes:
# (1) Effective drug combinations: Synergistic drug combinations from FimmDrugComb
#     tested on the cell lines for the particular disease.
# (2) Adverse drug combinations: DrugBank drug-drug interaction pairs with
#     risk/severity indications that contain both drugs involved
#     in forming effective pairs
# (3) Drug pairs also filtered to keep only those for which there are reported targets






# Load libraries
library(optparse)



# Get arguments
option_list = list(
  make_option(c("--disease"), type = "character", default = NULL, 
              help = "Name of the disease. The disease name will also be used as file name. e.g.: LungCancer, BreastCancer, etc.", metavar = "character")
)

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

if(is.null(opt$disease)){
  print_help(opt_parser)
  stop("--disease argument needed", call.=FALSE)
}

# Define global options for this script 
disease <- opt$disease





# Read drug-drug interactions from DrugBank involving risk or severity for side effects
DrugBank_drugInteractions <- readRDS("InputFiles/ReferenceList/DrugBank_drugInteractions_withRiskSeverity.rds")


# Read drug target interactions from drug bank
DrugBank_Drug_Target_Net <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")
Drug_Target_Net <- DrugBank_Drug_Target_Net





# Generate effective combinations ----------------------------------------------

# Read synergistic drug combinations from DrugCombDb for the specific disease
FimmDrugComb_drugCombs <- readRDS(paste0("InputFiles/ReferenceList/FimmDrugComb_", disease, "_drugCombinations.rds"))
effectiveCombinations <- FimmDrugComb_drugCombs$synergy 
print(paste0("Number of synergy combinations from DrugCombDb: ", nrow(effectiveCombinations)))


effectiveCombinations <- effectiveCombinations[effectiveCombinations$Drug1_DrugBank_drug_id != effectiveCombinations$Drug2_DrugBank_drug_id, ]
print(paste0("Number of synergy combinations after same drug removal: ", nrow(effectiveCombinations)))


# Remove drug pairs that are already reported as adverse pairs in DrugBank
remove_index <- c()
for(i in 1:nrow(effectiveCombinations)){
  drug1 <- effectiveCombinations[i, "Drug1_DrugBank_drug_id"] 
  drug2 <- effectiveCombinations[i, "Drug2_DrugBank_drug_id"] 
  
  tmp <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id == drug1 &  DrugBank_drugInteractions$Drug2_DrugBank_drug_id == drug2)|
                                     (DrugBank_drugInteractions$Drug1_DrugBank_drug_id == drug2 &  DrugBank_drugInteractions$Drug2_DrugBank_drug_id == drug1),] 
  if(nrow(tmp) > 0){remove_index <- c(remove_index, i)} 
} 
if(length(remove_index) > 0){effectiveCombinations <- effectiveCombinations[-remove_index,]}
rm(tmp)
print(paste0("Number of synergy combinations after adverse pair removal: ", nrow(effectiveCombinations)))


# Check if the drugs forming the pairs have reported targets
effectiveCombinations <- effectiveCombinations[(effectiveCombinations$Drug1_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id 
                                                & effectiveCombinations$Drug2_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id), ]
effectiveCombinations <- unique(effectiveCombinations)

print(paste0("Number of synergy combinations after target check: ", nrow(effectiveCombinations)))


effectiveCombinations_drugs <- unique(c(effectiveCombinations$Drug1_DrugBank_drug_id, effectiveCombinations$Drug2_DrugBank_drug_id))
print(paste0("Number of drugs forming effective combinations: ", length(effectiveCombinations_drugs)))





# Generate adverse combinations ------------------------------------------------

# Select DrugBank drug-drug interaction pairs that contain both drug 
# involved in forming effective pair
adverseCombinations <- DrugBank_drugInteractions[((DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in% effectiveCombinations_drugs) & 
                                                    (DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in% effectiveCombinations_drugs)),
                                                 c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")] 
adverseCombinations <- unique(adverseCombinations)


# Check if the drugs forming the pairs have reported targets
adverseCombinations <- adverseCombinations[(adverseCombinations$Drug1_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id 
                                            & adverseCombinations$Drug2_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id), ]





# Save drug combinations to file -----------------------------------------------
drugCombs <- list()
drugCombs$effectiveCombinations <- effectiveCombinations
drugCombs$adverseCombinations <- adverseCombinations


rownames(drugCombs$effectiveCombinations) <- rownames(drugCombs$adverseCombinations) <- NULL


if(!dir.exists("InputFiles/DrugCombinations/DrugCombs_v5/")){
  dir.create("InputFiles/DrugCombinations/DrugCombs_v5/", recursive = TRUE)
} 
saveRDS(drugCombs, paste0("InputFiles/DrugCombinations/DrugCombs_v5//DrugComb_", disease, "_v5.rds"))

cat(paste0("\n\n\nNumber of drug combinations for ", disease, ":\n"))
print(lapply(drugCombs, nrow))

print(warnings())