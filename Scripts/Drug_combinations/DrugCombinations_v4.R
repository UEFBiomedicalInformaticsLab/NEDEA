set.seed(5081)
rm(list = ls())



# Drug combinations (version 4)



# Notes:
# (1) Effective drug combinations: Selected drug pairs from
#     FimmDrugComb that were tested on disease related cell lines
#     and found synergistic
# (2) Adverse drug combinations: Selected drug pairs reported in DrugBank
#     to have risk or severity to adverse effect and containing both as approved drugs for the disease.
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





# Define the disease IDs
diseaseIds <- switch(disease,
                     "LungCancer" = c("EFO_0001071", "EFO_0003060", "EFO_0000702"),
                     "Leukemia" = c("EFO_0000565", "MONDO_0005059"),
                     "BreastCancer" = c("MONDO_0007254", "MONDO_0000618", "EFO_0005537", "EFO_0000305", "EFO_0000306"),
                     "ProstateCancer" = c("EFO_0000196", "MONDO_0008315", "EFO_0001663", "EFO_0000673", "EFO_1000499"),
                     "ColonCancer" = c("EFO_1001950", "EFO_1001949", "MONDO_0003978", "MONDO_0018513", "MONDO_0024479"),
                     "LiverCancer" = c("MONDO_0003243", "EFO_0000182", "MONDO_0002691", "MONDO_0004695", "MONDO_0002397"),
                     "OvaryCancer" = c("MONDO_0008170", "MONDO_0003812", "MONDO_0000548", "EFO_0001075", "EFO_1000416"),
                     "SkinCancer" = c("MONDO_0002898", "EFO_0009259", "MONDO_0002529", "EFO_1001471", "EFO_1000531"),
                     "KidneyCancer" = c("MONDO_0002367",  "EFO_0002890", "EFO_0003865"),
                     stop(paste("Disease IDs not defined for ", disease)))

cat(paste0("\n\nDisease IDs for ", disease, ": ", paste(diseaseIds, collapse = ", "), "\n\n"))



# Read drug-drug interactions from DrugBank involving risk or severity for side effects
DrugBank_drugInteractions <- readRDS("InputFiles/ReferenceList/DrugBank_drugInteractions_withRiskSeverity.rds")


# Read drug target interactions
DrugBank_Drug_Target_Net <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")
Drug_Target_Net <- DrugBank_Drug_Target_Net




# Retrieve all approved drugs for the disease
OpenTargets_Drug_Disease_Net <- readRDS("InputFiles/Associations/OpenTargets_Drug_Disease_Net.rds")
approved_drugs <- OpenTargets_Drug_Disease_Net[(OpenTargets_Drug_Disease_Net$Node2_disease_id %in% diseaseIds 
                                                & OpenTargets_Drug_Disease_Net$Drug_Disease_clinical_status == "Approved"), ]
approved_drugs <- unique(approved_drugs$Node1_drugbank_drug_id)
cat(paste0("\n\nNumber of approved drugs for ", disease, ": ", length(unique(approved_drugs)), "\n\n"))





# Generate adverse combinations ------------------------------------------------

# Select DrugBank drug-drug interaction pairs that contain at least one drug 
# involved in forming effective pair
adverseCombinations <- DrugBank_drugInteractions[((DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in% approved_drugs) & 
                                                    (DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in% approved_drugs)),
                                                 c("Drug1_DrugBank_drug_id", "Drug2_DrugBank_drug_id")] 


# Check if the drugs forming the pairs have reported targets
adverseCombinations <- adverseCombinations[(adverseCombinations$Drug1_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id 
                                            & adverseCombinations$Drug2_DrugBank_drug_id %in% Drug_Target_Net$Node1_drugbank_drug_id), ]





# Generate effective combinations ----------------------------------------------


# Read synergistic drug combinations from FimmDrugComb for the specific disease
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

print(paste0("Number of synergy combinations after target check: ", nrow(effectiveCombinations)))




# Save drug combinations to file -----------------------------------------------
drugCombs <- list()
drugCombs$effectiveCombinations <- effectiveCombinations
drugCombs$adverseCombinations <- adverseCombinations


rownames(drugCombs$effectiveCombinations) <- rownames(drugCombs$adverseCombinations) <- NULL


if(!dir.exists("InputFiles/DrugCombinations/DrugCombs_v4/")){
  dir.create("InputFiles/DrugCombinations/DrugCombs_v4/", recursive = TRUE)
} 
saveRDS(drugCombs, paste0("InputFiles/DrugCombinations/DrugCombs_v4//DrugComb_", disease, "_v4.rds"))

cat(paste0("\n\n\nNumber of drug combinations for ", disease, ":\n"))
print(lapply(drugCombs, nrow))

print(warnings())