<<<<<<< HEAD
set.seed(5081)



=======
>>>>>>> dd66bdd55a3da78129090252ed59959b311d68ad
# Enrichment analysis libraries for KidneyCancer


# Notes:
# Includes: renal carcinoma, kidney neoplasm, kidney cancer
# 
# NCBI MedGene Concept ID: C1378703, C0022665, C4682862 
#
# Experimental Factor Ontology: "EFO_0002890", "EFO_0003865"
# Mondo Disease Ontology: "MONDO_0002367"
#
# For Enrichr, used term ""





# Load libraries
library(org.Hs.eg.db)





# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_KidneyCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C1378703|C0022665|C4682862", 
                                                                 x = names(DisGeNET_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
names(DisGeNET_KidneyCancer2Gene_lib) <- paste0(names(DisGeNET_KidneyCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_KidneyCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_0002890|EFO_0003865|MONDO_0002367", 
                                                                             x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                             ignore.case = TRUE)]
names(OpenTargets_KidneyCancer2Gene_GA_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_GA_lib), "[OpenTargets_GA]")


# OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
# OpenTargets_KidneyCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_0002890|EFO_0003865|MONDO_0002367", 
#                                                                                x = names(OpenTargets_Disease2Gene_RNA_lib), 
#                                                                                ignore.case = TRUE)]
# names(OpenTargets_KidneyCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_KidneyCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_0002890|EFO_0003865|MONDO_0002367", 
                                                                               x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                               ignore.case = TRUE)]
names(OpenTargets_KidneyCancer2Gene_lit_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_lit_lib), "[OpenTargets_literature]")


# ## Enrichr
# 
# Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")
# 
# Enrichr_GeoKidneyCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "", 
#                                                                                        x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
#                                                                                        ignore.case = TRUE)]
# names(Enrichr_GeoKidneyCancerSignatures_Up) <- paste0(names(Enrichr_GeoKidneyCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")
# 
# Enrichr_GeoKidneyCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "", 
#                                                                                            x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
#                                                                                            ignore.case = TRUE)]
# names(Enrichr_GeoKidneyCancerSignatures_Down) <- paste0(names(Enrichr_GeoKidneyCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_KidneyCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Renal",
                                                               x = names(Intogen_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
names(Intogen_KidneyCancer2Gene_lib) <- paste0(names(Intogen_KidneyCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_KidneyCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Carcinoma, Renal Cell", 
                                                                 x = names(PharmGKB_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
names(PharmGKB_KidneyCancer2Gene_lib) <- paste0(names(PharmGKB_KidneyCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
KidneyCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0002890")]]


KidneyCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0002890"),]

KidneyCancer_Tscores <- tissue.specific.scores(disease_genes = KidneyCancer_genes$entrez, 
                                               ppi_network = ppi_strdb_700, 
                                               directed_network = FALSE, 
                                               tissue_expr_data = gtexv7_zscore,
                                               dis_relevant_tissues = KidneyCancer_rel_tissue_scores, 
                                               W = centrality_score$borda.disc, selected_tissues = NULL, 
                                               cutoff = 4, verbose = TRUE)

ThETA_KidneyCancer2Gene_lib <- row.names(KidneyCancer_Tscores[order(KidneyCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_KidneyCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_KidneyCancer2Gene_lib, c("ensembl_id")])
names(ThETA_KidneyCancer2Gene_lib) <- paste0("Renal carcinoma (EFO:0002890)", "[ThETA]")



# Merge and save into one

Enrichment_KidneyCancer2Gene_lib <- c(DisGeNET_KidneyCancer2Gene_lib, OpenTargets_KidneyCancer2Gene_GA_lib, 
                                      OpenTargets_KidneyCancer2Gene_lit_lib, 
                                      Intogen_KidneyCancer2Gene_lib,
                                      PharmGKB_KidneyCancer2Gene_lib, ThETA_KidneyCancer2Gene_lib)





# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_KidneyCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_KidneyCancer_lib.rds")



print(warnings())