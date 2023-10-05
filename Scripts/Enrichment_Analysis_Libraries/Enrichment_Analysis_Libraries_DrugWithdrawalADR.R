set.seed(5081)



# Enrichment libraries based on common ADRs related to drug withrawal
# 
# Primary terms include: hepatotoxicity, immune-related reaction, cardiotoxicity, 
# neurotoxicity, haematological toxicity, carcinogenicity 
# Ref: ADRs related to withdrawn drugs (Aronson, Jeffrey K. "Post-marketing drug withdrawals: 
# pharmacovigilance success, regulatory problems." Therapies 72.5 (2017): 555-561.)
# 
# 
# Additional terms manually extracted

# names(OpenTargets_Disease2Gene_lib)[grep("Neurotoxicity", names(OpenTargets_Disease2Gene_lib), ignore.case = TRUE)]

# Terms removed due to similarity
# - Cardiotoxicity (12.03.01.007) [ADReCS_gene]
# - Neurotoxicity (12.03.01.011) [ADReCS_gene]
# - Hepatitis, Drug-Induced (C1262760) [DisGeNET]
# - Drug-Induced Acute Liver Injury (C3658290) [DisGeNET]
# - Drug-Induced Liver Disease (C0860207) [DisGeNET]
# - Chemically-Induced Liver Toxicity (C4279912) [DisGeNET]





# Load libraries
library(openxlsx)



# Extract ADRs from ADReCS (ADR-gene associations) -----------------------------------------------------------

ADR_terms <- c("Adverse reaction (08.06.01.018)",
               "Cardiotoxicity (02.01.01.002)",
               "Neurotoxicity (17.02.10.002)",
               "Poisoning (12.03.01.004)",
               "Poisoning and toxicity (12.03.01)")
ADReCS_ADR2Gene_level3_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level3_lib.rds")
ADReCS_ADR2Gene_level4_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Gene_level4_lib.rds")
ADReCS_ADR2Gene_lib <- c(ADReCS_ADR2Gene_level4_lib, ADReCS_ADR2Gene_level3_lib)

ADReCS_gene_drugWithdrawalAdr_lib <- ADReCS_ADR2Gene_lib[names(ADReCS_ADR2Gene_lib) %in% ADR_terms]
names(ADReCS_gene_drugWithdrawalAdr_lib) <- paste0(names(ADReCS_gene_drugWithdrawalAdr_lib), " [ADReCS_gene]")



# Extract ADRs from ADReCS (ADR-protein associations) -----------------------------------------------------------

ADR_terms <- c("Adverse drug reaction (08.06.01.009)",
               "Cardiotoxicity (02.01.01.002)",
               "Drug toxicity (12.03.01.002)",
               "Haematotoxicity (01.05.01.007)",
               "Hepatotoxicity (09.01.07.009)",
               "Mitochondrial toxicity (12.03.01.009)",
               "Nephropathy toxic (12.03.01.010)",
               "Neurotoxicity (17.02.10.002)",
               "Ototoxicity (04.03.01.004)",
               "Poisoning and toxicity (12.03.01)")
ADReCS_ADR2Protein_level3_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Protein_level3_lib.rds")
ADReCS_ADR2Protein_level4_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/ADReCS_ADR2Protein_level4_lib.rds")
ADReCS_ADR2Protein_lib <- c(ADReCS_ADR2Protein_level4_lib, ADReCS_ADR2Protein_level3_lib)

ADReCS_protein_drugWithdrawalAdr_lib <- ADReCS_ADR2Protein_lib[names(ADReCS_ADR2Protein_lib) %in% ADR_terms]
names(ADReCS_protein_drugWithdrawalAdr_lib) <- paste0(names(ADReCS_protein_drugWithdrawalAdr_lib), " [ADReCS_protein]")



# Extract ADRs from CTD -----------------------------------------------------------

ADR_terms <- c("Cardiotoxicity (MESH:D066126)",
               "Chemical and Drug Induced Liver Injury (MESH:D056486)",
               "Chemical and Drug Induced Liver Injury, Chronic (MESH:D056487)",
               "Drug-Related Side Effects and Adverse Reactions (MESH:D064420)",
               "Neurotoxicity Syndromes (MESH:D020258)")
CTD_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/CTD_Disease2Gene_lib.rds")
CTD_drugWithdrawalAdr_lib <- CTD_Disease2Gene_lib[names(CTD_Disease2Gene_lib) %in% ADR_terms]
names(CTD_drugWithdrawalAdr_lib) <- paste0(names(CTD_drugWithdrawalAdr_lib), " [CTD]")



# Extract ADRs from DisGeNET -----------------------------------------------------------

ADR_terms <- c("Cardiotoxicity (C0876994)",
               "Chemical and Drug Induced Liver Injury (C4277682)",
               "Drug toxicity (C0013221)",
               "Neurotoxicity Syndromes (C0235032)")
DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_DiseaseGroup2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_DiseaseGroup2Gene_lib.rds")
DisGeNET_Phenotype2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Phenotype2Gene_lib.rds")
DisGeNET_Disease2Gene_lib <- c(DisGeNET_Disease2Gene_lib, DisGeNET_DiseaseGroup2Gene_lib, DisGeNET_Phenotype2Gene_lib)

DisGeNET_drugWithdrawalAdr_lib <- DisGeNET_Disease2Gene_lib[names(DisGeNET_Disease2Gene_lib) %in% ADR_terms]
names(DisGeNET_drugWithdrawalAdr_lib) <- paste0(names(DisGeNET_drugWithdrawalAdr_lib), " [DisGeNET]")



# Extract ADRs from OpenTargets -----------------------------------------------------------

ADR_terms <- c("ototoxicity (EFO_0006951)") # No relevant terms with good representation of genes
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_Disease2Gene_lib <- c(OpenTargets_Disease2Gene_GA_lib, OpenTargets_Disease2Gene_lit_lib, OpenTargets_Disease2Gene_RNA_lib)

OpenTargets_drugWithdrawalAdr_lib <- OpenTargets_Disease2Gene_lib[names(OpenTargets_Disease2Gene_lib) %in% ADR_terms]
names(OpenTargets_drugWithdrawalAdr_lib) <- paste0(names(OpenTargets_drugWithdrawalAdr_lib), " [OpenTargets]")



# Extract ADRs from PharmGKB -----------------------------------------------------------

ADR_terms <- c("cardiotoxicity",
               "Drug Toxicity",
               "drug-induced liver injury",
               "hematotoxicity",
               "Neurotoxicity Syndromes",
               "nephrotoxicity",
               "Ototoxicity")
PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")
PharmGKB_drugWithdrawalAdr_lib <- PharmGKB_Disease2Gene_lib[names(PharmGKB_Disease2Gene_lib) %in% ADR_terms]
names(PharmGKB_drugWithdrawalAdr_lib) <- paste0(names(PharmGKB_drugWithdrawalAdr_lib), " [PharmGKB]")



# Compile all gene list into one list
drugWithdrawal_Adr2Gene_lib <- c(ADReCS_gene_drugWithdrawalAdr_lib,
                                 ADReCS_protein_drugWithdrawalAdr_lib,
                                 CTD_drugWithdrawalAdr_lib,
                                 DisGeNET_drugWithdrawalAdr_lib,
                                 OpenTargets_drugWithdrawalAdr_lib,
                                 PharmGKB_drugWithdrawalAdr_lib)

if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(drugWithdrawal_Adr2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/drugWithdrawal_Adr2Gene_lib.rds")





# Check similarity between ADR terms
lib_size <- as.data.frame(lengths(drugWithdrawal_Adr2Gene_lib))
colnames(lib_size) <- "Number of genes"
lib_size$ADR <- row.names(lib_size)
lib_size <- lib_size[, c("ADR", "Number of genes")]

# Use Jaccard index to filter similar terms
lib_term_similarity <- data.frame()
for (i in unique(names(drugWithdrawal_Adr2Gene_lib))) {
  for (j in unique(names(drugWithdrawal_Adr2Gene_lib))) {
    if (i != j) {
      intersection <-
        intersect(drugWithdrawal_Adr2Gene_lib[[i]], drugWithdrawal_Adr2Gene_lib[[j]])
      union <- union(drugWithdrawal_Adr2Gene_lib[[i]], drugWithdrawal_Adr2Gene_lib[[j]])
      jaccard <- length(intersection) / length(union)
      lib_term_similarity <- rbind(lib_term_similarity,
                                   data.frame(
                                     ADR_1 = i,
                                     ADR_2 = j ,
                                     Jaccard = jaccard))
    }
  }
}
lib_term_similarity <- lib_term_similarity[order(lib_term_similarity$Jaccard, decreasing = TRUE), ]

if(!dir.exists("OutputFiles/Tables/")){ dir.create("OutputFiles/Tables/", recursive = TRUE) }
write.xlsx(list(Library_size = lib_size,
                Library_term_similarity = lib_term_similarity), 
           "OutputFiles/Tables/drugWithdrawal_Adr2Gene_libInfo.xlsx", overwrite = TRUE)



print(warnings())