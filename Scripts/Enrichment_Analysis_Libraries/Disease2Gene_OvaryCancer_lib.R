set.seed(5081)



# Enrichment analysis libraries for OvaryCancer


# Notes:
# Includes: ovarian carcinoma, Ovarian Endometrioid Adenocarcinoma, 
# ovarian cancer, ovarian endometrial cancer, ovarian clear cell cancer


# NCBI MedGene Concept ID: C2212006, C3544205, C0948216, C0346163, 
# C1335167, C1335177, C1140680, C0029925, C4721610
#
# Experimental Factor Ontology: , "EFO_0001075", "EFO_1000416"
# Mondo Disease Ontology: "MONDO_0008170", "MONDO_0003812", "MONDO_0000548"
# For Enrichr, used term ""





# Load libraries
library(org.Hs.eg.db)





# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_OvaryCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C2212006|C3544205|C0948216|C0346163|C1335167|C1335177|C1140680|C0029925|C4721610", 
                                                                 x = names(DisGeNET_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
names(DisGeNET_OvaryCancer2Gene_lib) <- paste0(names(DisGeNET_OvaryCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_OvaryCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_0001075|EFO_1000416|MONDO_0008170|MONDO_0003812|MONDO_0000548", 
                                                                             x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                             ignore.case = TRUE)]
names(OpenTargets_OvaryCancer2Gene_GA_lib) <- paste0(names(OpenTargets_OvaryCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_OvaryCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_0001075|EFO_1000416|MONDO_0008170|MONDO_0003812|MONDO_0000548", 
                                                                               x = names(OpenTargets_Disease2Gene_RNA_lib), 
                                                                               ignore.case = TRUE)]
names(OpenTargets_OvaryCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_OvaryCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_OvaryCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_0001075|EFO_1000416|MONDO_0008170|MONDO_0003812|MONDO_0000548", 
                                                                               x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                               ignore.case = TRUE)]
names(OpenTargets_OvaryCancer2Gene_lit_lib) <- paste0(names(OpenTargets_OvaryCancer2Gene_lit_lib), "[OpenTargets_literature]")


# ## Enrichr [NOT AVAILABLE]
# 
# Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")
# 
# Enrichr_GeoOvaryCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Breast Cancer", 
#                                                                                        x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
#                                                                                        ignore.case = TRUE)]
# names(Enrichr_GeoOvaryCancerSignatures_Up) <- paste0(names(Enrichr_GeoOvaryCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")
# 
# Enrichr_GeoOvaryCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Breast Cancer", 
#                                                                                            x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
#                                                                                            ignore.case = TRUE)]
# names(Enrichr_GeoOvaryCancerSignatures_Down) <- paste0(names(Enrichr_GeoOvaryCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_OvaryCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Ovary cancer",
                                                               x = names(Intogen_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
names(Intogen_OvaryCancer2Gene_lib) <- paste0(names(Intogen_OvaryCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_OvaryCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Ovarian Neoplasms", 
                                                                 x = names(PharmGKB_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
names(PharmGKB_OvaryCancer2Gene_lib) <- paste0(names(PharmGKB_OvaryCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
OvaryCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0001075")]]


OvaryCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0001075"),]

OvaryCancer_Tscores <- tissue.specific.scores(disease_genes = OvaryCancer_genes$entrez, 
                                               ppi_network = ppi_strdb_700, 
                                               directed_network = FALSE, 
                                               tissue_expr_data = gtexv7_zscore,
                                               dis_relevant_tissues = OvaryCancer_rel_tissue_scores, 
                                               W = centrality_score$borda.disc, selected_tissues = NULL, 
                                               cutoff = 4, verbose = TRUE)

ThETA_OvaryCancer2Gene_lib <- row.names(OvaryCancer_Tscores[order(OvaryCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_OvaryCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_OvaryCancer2Gene_lib, c("ensembl_id")])
names(ThETA_OvaryCancer2Gene_lib) <- paste0("Ovarian carcinoma (EFO_0001075)", "[ThETA]")



# Merge and save into one

Enrichment_OvaryCancer2Gene_lib <- c(DisGeNET_OvaryCancer2Gene_lib, OpenTargets_OvaryCancer2Gene_GA_lib, 
                                      OpenTargets_OvaryCancer2Gene_RNA_lib, OpenTargets_OvaryCancer2Gene_lit_lib, 
                                      Intogen_OvaryCancer2Gene_lib,
                                      PharmGKB_OvaryCancer2Gene_lib, ThETA_OvaryCancer2Gene_lib)





# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_OvaryCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_OvaryCancer_lib.rds")



print(warnings())