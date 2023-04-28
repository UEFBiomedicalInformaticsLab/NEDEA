# Enrichment analysis libraries for SkinCancer


# Notes:
# Includes: skin carcinoma, Merkel cell skin cancer, Skin Sarcoma, skin cancer, skin squamous cell carcinoma

# NCBI MedGene Concept ID: C1321872|C0007114 

# Experimental Factor Ontology:  "EFO_0009259", "EFO_1001471", "EFO_1000531", 
# Mondo Disease Ontology: "MONDO_0002898", "MONDO_0002529"
# For Enrichr, used term 





# Load libraries
library(org.Hs.eg.db)





# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_SkinCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep(pattern = "C1321872|C0007114", 
                                                                x = names(DisGeNET_Disease2Gene_lib), 
                                                                ignore.case = TRUE)]
names(DisGeNET_SkinCancer2Gene_lib) <- paste0(names(DisGeNET_SkinCancer2Gene_lib), "[DisGeNET_curated]")


## OpenTargets
OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_SkinCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep(pattern = "EFO_0009259|EFO_1001471|EFO_1000531|MONDO_0002898|MONDO_0002529", 
                                                                            x = names(OpenTargets_Disease2Gene_GA_lib), 
                                                                            ignore.case = TRUE)]
names(OpenTargets_SkinCancer2Gene_GA_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_GA_lib), "[OpenTargets_GA]")


# OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
# OpenTargets_SkinCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep(pattern = "EFO_0009259|EFO_1001471|EFO_1000531|MONDO_0002898|MONDO_0002529",
#                                                                                  x = names(OpenTargets_Disease2Gene_RNA_lib),
#                                                                                  ignore.case = TRUE)]
# names(OpenTargets_SkinCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_SkinCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep(pattern = "EFO_0009259|EFO_1001471|EFO_1000531|MONDO_0002898|MONDO_0002529", 
                                                                              x = names(OpenTargets_Disease2Gene_lit_lib), 
                                                                              ignore.case = TRUE)]
names(OpenTargets_SkinCancer2Gene_lit_lib) <- paste0(names(OpenTargets_SkinCancer2Gene_lit_lib), "[OpenTargets_literature]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoSkinCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep(pattern = "In situ melanoma of skin", 
                                                                                      x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                                                      ignore.case = TRUE)]
names(Enrichr_GeoSkinCancerSignatures_Up) <- paste0(names(Enrichr_GeoSkinCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoSkinCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep(pattern = "In situ melanoma of skin", 
                                                                                          x = names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                                          ignore.case = TRUE)]
names(Enrichr_GeoSkinCancerSignatures_Down) <- paste0(names(Enrichr_GeoSkinCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")


## Intogen
Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/Intogen_Disease2Gene_lib.rds")
Intogen_SkinCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep(pattern = "skin",
                                                              x = names(Intogen_Disease2Gene_lib), 
                                                              ignore.case = TRUE)]
names(Intogen_SkinCancer2Gene_lib) <- paste0(names(Intogen_SkinCancer2Gene_lib), "[Intogen]")


# ## PharmGKB [NONE FOUND]
# 
# PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_Analysis_Libraries/PharmGKB_Disease2Gene_lib.rds")
# 
# PharmGKB_SkinCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep(pattern = "", 
#                                                                 x = names(PharmGKB_Disease2Gene_lib), 
#                                                                 ignore.case = TRUE)]
# names(PharmGKB_SkinCancer2Gene_lib) <- paste0(names(PharmGKB_SkinCancer2Gene_lib), "[PharmGKB_associated]")


# ## ThETA [NONE FOUND]
# library(ThETA)
# source("Scripts/ThETA/Corrected_Functions.R")
# data(gtexv7_zscore)
# data(ppi_strdb_700)
# data(dis_vrnts)
# data(disease_tissue_zscores)
# data(centrality_score)
# 
# #  compile the tissue-specific efficacy estimates of target(gene)-disease associations
# SkinCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "")]]
# 
# 
# SkinCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == ""),]
# 
# SkinCancer_Tscores <- tissue.specific.scores(disease_genes = SkinCancer_genes$entrez, 
#                                               ppi_network = ppi_strdb_700, 
#                                               directed_network = FALSE, 
#                                               tissue_expr_data = gtexv7_zscore,
#                                               dis_relevant_tissues = SkinCancer_rel_tissue_scores, 
#                                               W = centrality_score$borda.disc, selected_tissues = NULL, 
#                                               cutoff = 4, verbose = TRUE)
# 
# ThETA_SkinCancer2Gene_lib <- row.names(SkinCancer_Tscores[order(SkinCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
# ThETA_SkinCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_SkinCancer2Gene_lib, c("ensembl_id")])
# names(ThETA_SkinCancer2Gene_lib) <- paste0("", "[ThETA]")



# Merge and save into one

Enrichment_SkinCancer2Gene_lib <- c(DisGeNET_SkinCancer2Gene_lib, OpenTargets_SkinCancer2Gene_GA_lib, 
                                     # OpenTargets_SkinCancer2Gene_RNA_lib, 
                                     OpenTargets_SkinCancer2Gene_lit_lib, 
                                     Enrichr_GeoSkinCancerSignatures_Up, Enrichr_GeoSkinCancerSignatures_Down,
                                     Intogen_SkinCancer2Gene_lib)
                                     # PharmGKB_SkinCancer2Gene_lib, 
                                     # ThETA_SkinCancer2Gene_lib)




# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(Enrichment_SkinCancer2Gene_lib, "InputFiles/Enrichment_Analysis_Libraries/Disease2Gene_SkinCancer_lib.rds")



print(warnings())