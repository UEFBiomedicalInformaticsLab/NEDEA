set.seed(5081)



# Enrichment analysis libraries for ProstateCancer



# Load libraries
library(org.Hs.eg.db)



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")



## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_ProstateCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep("prostate", 
                                         names(DisGeNET_Disease2Gene_lib), 
                                         ignore.case = TRUE)]
DisGeNET_ProstateCancer2Gene_lib <- DisGeNET_ProstateCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                      names(DisGeNET_ProstateCancer2Gene_lib), 
                      ignore.case = TRUE)]
DisGeNET_ProstateCancer2Gene_lib <- DisGeNET_ProstateCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                      names(DisGeNET_ProstateCancer2Gene_lib), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(DisGeNET_ProstateCancer2Gene_lib) <- paste0(names(DisGeNET_ProstateCancer2Gene_lib), "[DisGeNET_curated]")



## OpenTargets

OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_ProstateCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep("prostate", 
                                               names(OpenTargets_Disease2Gene_GA_lib), 
                                               ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_GA_lib <- OpenTargets_ProstateCancer2Gene_GA_lib[grep("cancer|carcinoma|sarcoma", 
                      names(OpenTargets_ProstateCancer2Gene_GA_lib), 
                      ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_GA_lib <- OpenTargets_ProstateCancer2Gene_GA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                      names(OpenTargets_ProstateCancer2Gene_GA_lib), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(OpenTargets_ProstateCancer2Gene_GA_lib) <- paste0(names(OpenTargets_ProstateCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_ProstateCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep("prostate", 
                                                names(OpenTargets_Disease2Gene_RNA_lib), 
                                                ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_RNA_lib <- OpenTargets_ProstateCancer2Gene_RNA_lib[grep("cancer|carcinoma|sarcoma", 
                      names(OpenTargets_ProstateCancer2Gene_RNA_lib), 
                      ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_RNA_lib <- OpenTargets_ProstateCancer2Gene_RNA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                      names(OpenTargets_ProstateCancer2Gene_RNA_lib), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(OpenTargets_ProstateCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_ProstateCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_ProstateCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep("prostate", 
                                                names(OpenTargets_Disease2Gene_lit_lib), 
                                                ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_lit_lib <- OpenTargets_ProstateCancer2Gene_lit_lib[grep("cancer|carcinoma|sarcoma", 
                      names(OpenTargets_ProstateCancer2Gene_lit_lib), 
                      ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_lit_lib <- OpenTargets_ProstateCancer2Gene_lit_lib[grep("hereditary|familial|susceptibility|predisposition", 
                      names(OpenTargets_ProstateCancer2Gene_lit_lib), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(OpenTargets_ProstateCancer2Gene_lit_lib) <- paste0(names(OpenTargets_ProstateCancer2Gene_lit_lib), "[OpenTargets_literature]")


OpenTargets_Disease2Gene_SM_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_SM_lib.rds")
OpenTargets_ProstateCancer2Gene_SM_lib <- OpenTargets_Disease2Gene_SM_lib[grep("prostate", 
                                                                               names(OpenTargets_Disease2Gene_SM_lib), 
                                                                               ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_SM_lib <- OpenTargets_ProstateCancer2Gene_SM_lib[grep("cancer|carcinoma|sarcoma", 
                                                                                      names(OpenTargets_ProstateCancer2Gene_SM_lib), 
                                                                                      ignore.case = TRUE)]
OpenTargets_ProstateCancer2Gene_SM_lib <- OpenTargets_ProstateCancer2Gene_SM_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                      names(OpenTargets_ProstateCancer2Gene_SM_lib), 
                                                                                      ignore.case = TRUE, 
                                                                                      invert = TRUE)] 
names(OpenTargets_ProstateCancer2Gene_SM_lib) <- paste0(names(OpenTargets_ProstateCancer2Gene_SM_lib), "[OpenTargets_SM]")


## Enrichr

Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")

Enrichr_GeoProstateCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("prostate", 
                                                            names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up), 
                                                            ignore.case = TRUE)]
Enrichr_GeoProstateCancerSignatures_Up <- Enrichr_GeoProstateCancerSignatures_Up[grep("cancer|carcinoma|sarcoma", 
                            names(Enrichr_GeoProstateCancerSignatures_Up), 
                            ignore.case = TRUE)]
Enrichr_GeoProstateCancerSignatures_Up <- Enrichr_GeoProstateCancerSignatures_Up[grep("hereditary|familial|susceptibility|predisposition", 
                      names(Enrichr_GeoProstateCancerSignatures_Up), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(Enrichr_GeoProstateCancerSignatures_Up) <- paste0(names(Enrichr_GeoProstateCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")

Enrichr_GeoProstateCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("prostate", 
                                                                names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down), 
                                                                ignore.case = TRUE)]
Enrichr_GeoProstateCancerSignatures_Down <- Enrichr_GeoProstateCancerSignatures_Down[grep("cancer|carcinoma|sarcoma", 
                                names(Enrichr_GeoProstateCancerSignatures_Down), 
                                ignore.case = TRUE)]
Enrichr_GeoProstateCancerSignatures_Down <- Enrichr_GeoProstateCancerSignatures_Down[grep("hereditary|familial|susceptibility|predisposition", 
                      names(Enrichr_GeoProstateCancerSignatures_Down), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(Enrichr_GeoProstateCancerSignatures_Down) <- paste0(names(Enrichr_GeoProstateCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")



## Intogen

Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Intogen_Disease2Gene_lib.rds")
Intogen_ProstateCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep("prostate", names(Intogen_Disease2Gene_lib), 
                                        ignore.case = TRUE)]
Intogen_ProstateCancer2Gene_lib <- Intogen_ProstateCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                      names(Intogen_ProstateCancer2Gene_lib), 
                      ignore.case = TRUE)]
Intogen_ProstateCancer2Gene_lib <- Intogen_ProstateCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                      names(Intogen_ProstateCancer2Gene_lib), 
                      ignore.case = TRUE, 
                      invert = TRUE)] 
names(Intogen_ProstateCancer2Gene_lib) <- paste0(names(Intogen_ProstateCancer2Gene_lib), "[Intogen]")



# ## PharmGKB
# 
# PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/PharmGKB_Disease2Gene_lib.rds")
# PharmGKB_ProstateCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep("prostate",
#                                                                    names(PharmGKB_Disease2Gene_lib),
#                                                                    ignore.case = TRUE)]
# PharmGKB_ProstateCancer2Gene_lib <- PharmGKB_ProstateCancer2Gene_lib[grep("cancer|carcinoma|sarcoma",
#                                                                           names(PharmGKB_ProstateCancer2Gene_lib),
#                                                                           ignore.case = TRUE)]
# PharmGKB_ProstateCancer2Gene_lib <- PharmGKB_ProstateCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition",
#                                                                           names(PharmGKB_ProstateCancer2Gene_lib),
#                                                                           ignore.case = TRUE,
#                                                                           invert = TRUE)]
# names(PharmGKB_ProstateCancer2Gene_lib) <- paste0(names(PharmGKB_ProstateCancer2Gene_lib), "[PharmGKB_associated]")



## ThETA
library(ThETA)
source("Scripts/ThETA/Corrected_Functions.R")
data(gtexv7_zscore)
data(ppi_strdb_700)
data(dis_vrnts)
data(disease_tissue_zscores)
data(centrality_score)

#  compile the tissue-specific efficacy estimates of target(gene)-disease associations
ProstateCancer_genes = dis_vrnts[[which(names(dis_vrnts) == "EFO:0001663")]]


ProstateCancer_rel_tissue_scores = disease_tissue_zscores$z[which(rownames(disease_tissue_zscores$z) == "EFO:0001663"),]

ProstateCancer_Tscores <- tissue.specific.scores(disease_genes = ProstateCancer_genes$entrez, 
                                                 ppi_network = ppi_strdb_700, 
                                                 directed_network = FALSE, 
                                                 tissue_expr_data = gtexv7_zscore,
                                                 dis_relevant_tissues = ProstateCancer_rel_tissue_scores, 
                                                 W = centrality_score$borda.disc, selected_tissues = NULL, 
                                                 cutoff = 4, verbose = TRUE)

ThETA_ProstateCancer2Gene_lib <- row.names(ProstateCancer_Tscores[order(ProstateCancer_Tscores$avg_tissue_score, decreasing = TRUE)[1:50],])
ThETA_ProstateCancer2Gene_lib <- list(entrezId_2_ensemblId[entrezId_2_ensemblId$gene_id %in% ThETA_ProstateCancer2Gene_lib, c("ensembl_id")])
names(ThETA_ProstateCancer2Gene_lib) <- paste0("Prostate carcinoma (EFO_0001663)", "[ThETA]")



# Merge and save into one

Enrichment_ProstateCancer2Gene_lib <- c(DisGeNET_ProstateCancer2Gene_lib, OpenTargets_ProstateCancer2Gene_GA_lib, 
                                        OpenTargets_ProstateCancer2Gene_RNA_lib, OpenTargets_ProstateCancer2Gene_lit_lib, 
                                        OpenTargets_ProstateCancer2Gene_SM_lib,
                                        Enrichr_GeoProstateCancerSignatures_Up, Enrichr_GeoProstateCancerSignatures_Down,
                                        Intogen_ProstateCancer2Gene_lib, ThETA_ProstateCancer2Gene_lib)



# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Enrichment_ProstateCancer2Gene_lib, "InputFiles/Enrichment_analysis_libraries/Disease2Gene_ProstateCancer_lib.rds")



print(warnings())