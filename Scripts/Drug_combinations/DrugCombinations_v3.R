set.seed(5081)
rm(list = ls())



# Drug combinations (version 3)



# Notes:
# (1) Used all the approved drugs for the disease to generate all possible combinations.
# (2) The pairs already reported to have adverse drug-drug interactions 
#     in DrugBank were considered adverse pairs else effective pairs






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
                     stop(paste("Disease IDs not defined for ", disease)))

cat(paste0("\n\nDisease IDs for ", disease, ": ", paste(diseaseIds, collapse = ", "), "\n\n"))



# Read drug-drug interactions from DrugBank involving risk or severity for side effects
DrugBank_drugInteractions <- readRDS("InputFiles/ReferenceList/DrugBank_drugInteractions_withRiskSeverity.rds")


# Read drug target interactions from drug bank
DrugBank_Drug_Target_Net <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_Net.rds")





# Retrieve all approved drugs for the disease
OpenTargets_Drug_Disease_Net <- readRDS("InputFiles/Associations/OpenTargets_Drug_Disease_Net.rds")
approved_drugs <- OpenTargets_Drug_Disease_Net[(OpenTargets_Drug_Disease_Net$Node2_disease_id %in% diseaseIds 
                                                & OpenTargets_Drug_Disease_Net$Drug_Disease_clinical_status == "Approved"), ]
approved_drugs <- unique(approved_drugs$Node1_drugbank_drug_id)
cat(paste0("\n\nNumber of approved drugs for ", disease, ": ", length(unique(approved_drugs)), "\n\n"))





# Generate the drug combinations
effectiveCombinations <- data.frame(Drug1_DrugBank_drug_id = character(), Drug2_DrugBank_drug_id = character())
adverseCombinations <- data.frame(Drug1_DrugBank_drug_id = character(), Drug2_DrugBank_drug_id = character())

for(i in 1:length(approved_drugs)){
  for(j in 1:length(approved_drugs)){
    
    if(j > i){
      
      drug1 <- approved_drugs[i]
      drug2 <- approved_drugs[j]
      
      
      tmp1 <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in%  drug1 &
                                           DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in%  drug2),  ] 
      tmp2 <- DrugBank_drugInteractions[(DrugBank_drugInteractions$Drug1_DrugBank_drug_id %in%  drug2 &
                                           DrugBank_drugInteractions$Drug2_DrugBank_drug_id %in%  drug1),  ] 
      
      tmp3 <- data.frame("Drug1_DrugBank_drug_id" = drug1, 
                         "Drug2_DrugBank_drug_id" = drug2)
      
      
      if(nrow(rbind(tmp1, tmp2))>0){
        adverseCombinations <- rbind(adverseCombinations, tmp3)
      }else{ effectiveCombinations <- rbind(effectiveCombinations, tmp3) } 
      
    }
  }
}


# Check if the drugs forming the pairs have reported targets
effectiveCombinations <- effectiveCombinations[(effectiveCombinations$Drug1_DrugBank_drug_id %in% DrugBank_Drug_Target_Net$Node1_drugbank_drug_id 
                                                & effectiveCombinations$Drug2_DrugBank_drug_id %in% DrugBank_Drug_Target_Net$Node1_drugbank_drug_id), ]

adverseCombinations <- adverseCombinations[(adverseCombinations$Drug1_DrugBank_drug_id %in% DrugBank_Drug_Target_Net$Node1_drugbank_drug_id 
                                            & adverseCombinations$Drug2_DrugBank_drug_id %in% DrugBank_Drug_Target_Net$Node1_drugbank_drug_id), ]





# Save drug combinations to file -----------------------------------------------
drugCombs <- list()
drugCombs$effectiveCombinations <- effectiveCombinations
drugCombs$adverseCombinations <- adverseCombinations


rownames(drugCombs$effectiveCombinations) <- rownames(drugCombs$adverseCombinations) <- NULL


if(!dir.exists("InputFiles/DrugCombinations/DrugCombs_v3/")){
  dir.create("InputFiles/DrugCombinations/DrugCombs_v3/", recursive = TRUE)
} 
saveRDS(drugCombs, paste0("InputFiles/DrugCombinations/DrugCombs_v3//DrugComb_", disease, "_v3.rds"))

cat(paste0("\n\n\nNumber of drug combinations for ", disease, ":\n"))
print(lapply(drugCombs, nrow))

print(warnings())