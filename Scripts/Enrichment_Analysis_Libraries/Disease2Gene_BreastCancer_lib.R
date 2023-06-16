<<<<<<< HEAD
set.seed(5081)



=======
>>>>>>> dd66bdd55a3da78129090252ed59959b311d68ad
# Enrichment analysis libraries for BreastCancer


# Notes:
# Includes: triple-negative breast cancer, breast carcinoma, 
# breast tumor luminal, breast cancer, Her2-receptor negative breast cancer

# NCBI MedGene Concept ID: C0678222, C0006142, C0858252, 
# C0278601, C1134719, C3539878, C1527349, C1336076
#
# Experimental Factor Ontology: "EFO_0005537", "EFO_0000305", "EFO_0000306"
# Mondo Disease Ontology: "MONDO_0007254", "MONDO_0000618"
# For Enrichr, used term ""





# Load libraries
library(org.Hs.eg.db)





# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_BreastCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C0678222|C0006142|C0858252|C0278601|C1134719|C3539878|C1527349|C1336076", 
                                                                x = names(DisGeNET_Disease2Gene_lib), 
                                                                ignore.case = TRUE)]
names(DisGeNET_BreastCancer2Gene_lib) <- paste0(names(DisGeNET_BreastCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_BreastCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "MONDO_0007254|MONDO_0000618|EFO_0005537|EFO_0000305|EFO_0000306", 
                                                                            x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                            ignore.case = TRUE)]
names(OpenTargets_BreastCancer2Gene_GA_lib) <- paste0(names(OpenTargets_BreastCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_BreastCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "MONDO_0007254|MONDO_0000618|EFO_0005537|EFO_0000305|EFO_0000306", 
                                                                              x = names(OpenTargets_Disease2Gene_RNA_lib), 
                                                                              ignore.case = TRUE)]
names(OpenTargets_BreastCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_BreastCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_BreastCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "MONDO_0007254|MONDO_0000618|EFO_0005537|EFO_0000305|EFO_0000306", 
                                                                              x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                              ignore.case = TRUE)]
names(OpenTargets_BreastCancer2Gene_lit_lib) <- paste0(names(OpenTargets_BreastCancer2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoBreastCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Breast Cancer", 
                                                                                      x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                      ignore.case = TRUE)]
names(Enrichr_GeoBreastCancerSignatures_Up) <- paste0(names(Enrichr_GeoBreastCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoBreastCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Breast Cancer", 
                                                                                          x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                          ignore.case = TRUE)]
names(Enrichr_GeoBreastCancerSignatures_Down) <- paste0(names(Enrichr_GeoBreastCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_BreastCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Breast adenocarcinoma",
                                                              x = names(Intogen_Disease2Gene_lib), 
                                                              ignore.case = TRUE)]
names(Intogen_BreastCancer2Gene_lib) <- paste0(names(Intogen_BreastCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_BreastCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Breast Neoplasms", 
                                                                x = names(PharmGKB_Disease2Gene_lib), 
                                                                ignore.case = TRUE)]
names(PharmGKB_BreastCancer2Gene_lib) <- paste0(names(PharmGKB_BreastCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
BreastCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0000305")]]


BreastCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0000305"),]

BreastCancer_Tscores <- tissue.specific.scores(disease_genes = BreastCancer_genes$entrez, 
                                              ppi_network = ppi_strdb_700, 
                                              directed_network = FALSE, 
                                              tissue_expr_data = gtexv7_zscore,
                                              dis_relevant_tissues = BreastCancer_rel_tissue_scores, 
                                              W = centrality_score$borda.disc, selected_tissues = NULL, 
                                              cutoff = 4, verbose = TRUE)

ThETA_BreastCancer2Gene_lib <- row.names(BreastCancer_Tscores[order(BreastCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_BreastCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_BreastCancer2Gene_lib, c("ensembl_id")])
names(ThETA_BreastCancer2Gene_lib) <- paste0("Breast carcinoma (EFO_0000305)", "[ThETA]")



# Merge and save into one

Enrichment_BreastCancer2Gene_lib <- c(DisGeNET_BreastCancer2Gene_lib, OpenTargets_BreastCancer2Gene_GA_lib, 
                                     OpenTargets_BreastCancer2Gene_RNA_lib, OpenTargets_BreastCancer2Gene_lit_lib, 
                                     Enrichr_GeoBreastCancerSignatures_Up, Enrichr_GeoBreastCancerSignatures_Down,
                                     Intogen_BreastCancer2Gene_lib,
                                     PharmGKB_BreastCancer2Gene_lib, ThETA_BreastCancer2Gene_lib)





# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_BreastCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_BreastCancer_lib.rds")

print(warnings())