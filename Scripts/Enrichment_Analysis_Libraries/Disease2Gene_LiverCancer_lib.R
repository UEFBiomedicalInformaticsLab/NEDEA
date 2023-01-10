set.seed(5081)
rm(list = ls())



# Enrichment analysis (Disease2Gene) library for Liver Cancer


# Notes:
# Includes: hepatocellular carcinoma, hepatocellular clear cell carcinoma, liver cancer,
# liver lymphoma, liver sarcoma

# NCBI MedGene Concept ID: C2239176
# Experimental Factor Ontology: "EFO_0000182"
# Mondo Disease Ontology: "MONDO_0003243", "MONDO_0002691", "MONDO_0004695", "MONDO_0002397"
# For Enrichr, used term "Carcinoma, Hepatocellular"




# Load libraries
library(org.Hs.eg.db)

# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_LiverCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C2239176", 
                                                                x = names(DisGeNET_Disease2Gene_lib), 
                                                                ignore.case = TRUE)]
names(DisGeNET_LiverCancer2Gene_lib) <- paste0(names(DisGeNET_LiverCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_LiverCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "MONDO_0003243|MONDO_0002691|MONDO_0004695|MONDO_0002397|EFO_0000182", 
                                                                            x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                            ignore.case = TRUE)]
names(OpenTargets_LiverCancer2Gene_GA_lib) <- paste0(names(OpenTargets_LiverCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_LiverCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "MONDO_0003243|MONDO_0002691|MONDO_0004695|MONDO_0002397|EFO_0000182", 
                                                                              x = names(OpenTargets_Disease2Gene_RNA_lib), 
                                                                              ignore.case = TRUE)]
names(OpenTargets_LiverCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_LiverCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_LiverCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "MONDO_0003243|MONDO_0002691|MONDO_0004695|MONDO_0002397|EFO_0000182", 
                                                                              x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                              ignore.case = TRUE)]
names(OpenTargets_LiverCancer2Gene_lit_lib) <- paste0(names(OpenTargets_LiverCancer2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoLiverCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "Carcinoma, Hepatocellular", 
                                                                                     x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                     ignore.case = TRUE)]
names(Enrichr_GeoLiverCancerSignatures_Up) <- paste0(names(Enrichr_GeoLiverCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoLiverCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "Carcinoma, Hepatocellular", 
                                                                                         x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                         ignore.case = TRUE)]
names(Enrichr_GeoLiverCancerSignatures_Down) <- paste0(names(Enrichr_GeoLiverCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_LiverCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "Hepatic cancer",
                                                             x = names(Intogen_Disease2Gene_lib), 
                                                             ignore.case = TRUE)]
names(Intogen_LiverCancer2Gene_lib) <- paste0(names(Intogen_LiverCancer2Gene_lib), "[Intogen]")


## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")

PharmGKB_LiverCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "Carcinoma, Hepatocellular", 
                                                               x = names(PharmGKB_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
names(PharmGKB_LiverCancer2Gene_lib) <- paste0(names(PharmGKB_LiverCancer2Gene_lib), "[PharmGKB_associated]")


## ThETA
library(ThETA)
source("ExternalTools/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  Compile the tissue-specific efficacy estimates of target(gene)-disease associations
LiverCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0000182")]]
LiverCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0000182"),]
LiverCancer_Tscores <- tissue.specific.scores(disease_genes = LiverCancer_genes$entrez, 
                                             ppi_network = ppi_strdb_700, 
                                             directed_network = FALSE, 
                                             tissue_expr_data = gtexv7_zscore,
                                             dis_relevant_tissues = LiverCancer_rel_tissue_scores, 
                                             W = centrality_score$borda.disc, selected_tissues = NULL, 
                                             cutoff = 4, verbose = TRUE)

ThETA_LiverCancer2Gene_lib <- row.names(LiverCancer_Tscores[order(LiverCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_LiverCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_LiverCancer2Gene_lib, c("ensembl_id")])
names(ThETA_LiverCancer2Gene_lib) <- paste0("hepatocellular carcinoma (EFO_0000182)", "[ThETA]")



# Merge and save into one

Enrichment_LiverCancer2Gene_lib <- c(DisGeNET_LiverCancer2Gene_lib, OpenTargets_LiverCancer2Gene_GA_lib, 
                                    OpenTargets_LiverCancer2Gene_RNA_lib, OpenTargets_LiverCancer2Gene_lit_lib, 
                                    Enrichr_GeoLiverCancerSignatures_Up, Enrichr_GeoLiverCancerSignatures_Down,
                                    Intogen_LiverCancer2Gene_lib,
                                    PharmGKB_LiverCancer2Gene_lib, ThETA_LiverCancer2Gene_lib)


# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_LiverCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_LiverCancer_lib.rds")

print(warnings())