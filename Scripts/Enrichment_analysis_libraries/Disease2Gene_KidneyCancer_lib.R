set.seed(5081)



# Enrichment analysis libraries for KidneyCancer



# Load libraries
library(org.Hs.eg.db)



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")



## DisGeNET

DisGeNET_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/DisGeNET_Disease2Gene_lib.rds")
DisGeNET_KidneyCancer2Gene_lib <- DisGeNET_Disease2Gene_lib[grep("kidney|renal", 
                                                                 names(DisGeNET_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
DisGeNET_KidneyCancer2Gene_lib <- DisGeNET_KidneyCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                      names(DisGeNET_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE)]
DisGeNET_KidneyCancer2Gene_lib <- DisGeNET_KidneyCancer2Gene_lib[grep("adrenal", 
                                                                      names(DisGeNET_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE,
                                                                      invert = TRUE)]
DisGeNET_KidneyCancer2Gene_lib <- DisGeNET_KidneyCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                      names(DisGeNET_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE, 
                                                                      invert = TRUE)] 
names(DisGeNET_KidneyCancer2Gene_lib) <- paste0(names(DisGeNET_KidneyCancer2Gene_lib), "[DisGeNET_curated]")



## OpenTargets

OpenTargets_Disease2Gene_GA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_GA_lib.rds")
OpenTargets_KidneyCancer2Gene_GA_lib <- OpenTargets_Disease2Gene_GA_lib[grep("kidney|renal", 
                                                                             names(OpenTargets_Disease2Gene_GA_lib), 
                                                                             ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_GA_lib <- OpenTargets_KidneyCancer2Gene_GA_lib[grep("cancer|carcinoma|sarcoma", 
                                                                                  names(OpenTargets_KidneyCancer2Gene_GA_lib), 
                                                                                  ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_GA_lib <- OpenTargets_KidneyCancer2Gene_GA_lib[grep("adrenal", 
                                                                                  names(OpenTargets_KidneyCancer2Gene_GA_lib), 
                                                                                  ignore.case = TRUE,
                                                                                  invert = TRUE)]
OpenTargets_KidneyCancer2Gene_GA_lib <- OpenTargets_KidneyCancer2Gene_GA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                  names(OpenTargets_KidneyCancer2Gene_GA_lib), 
                                                                                  ignore.case = TRUE, 
                                                                                  invert = TRUE)] 
names(OpenTargets_KidneyCancer2Gene_GA_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_GA_lib), "[OpenTargets_GA]")


OpenTargets_Disease2Gene_RNA_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_RNA_lib.rds")
OpenTargets_KidneyCancer2Gene_RNA_lib <- OpenTargets_Disease2Gene_RNA_lib[grep("kidney|renal",
                                                                               names(OpenTargets_Disease2Gene_RNA_lib),
                                                                               ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_RNA_lib <- OpenTargets_KidneyCancer2Gene_RNA_lib[grep("cancer|carcinoma|sarcoma",
                                                                                    names(OpenTargets_KidneyCancer2Gene_RNA_lib),
                                                                                    ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_RNA_lib <- OpenTargets_KidneyCancer2Gene_RNA_lib[grep("adrenal", 
                                                                                    names(OpenTargets_KidneyCancer2Gene_RNA_lib), 
                                                                                    ignore.case = TRUE,
                                                                                    invert = TRUE)]
OpenTargets_KidneyCancer2Gene_RNA_lib <- OpenTargets_KidneyCancer2Gene_RNA_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                    names(OpenTargets_KidneyCancer2Gene_RNA_lib), 
                                                                                    ignore.case = TRUE, 
                                                                                    invert = TRUE)] 
names(OpenTargets_KidneyCancer2Gene_RNA_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_RNA_lib), "[OpenTargets_RNA]")


OpenTargets_Disease2Gene_lit_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/OpenTargets_Disease2Gene_lit_lib.rds")
OpenTargets_KidneyCancer2Gene_lit_lib <- OpenTargets_Disease2Gene_lit_lib[grep("kidney|renal", 
                                                                               names(OpenTargets_Disease2Gene_lit_lib), 
                                                                               ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_lit_lib <- OpenTargets_KidneyCancer2Gene_lit_lib[grep("cancer|carcinoma|sarcoma", 
                                                                                    names(OpenTargets_KidneyCancer2Gene_lit_lib), 
                                                                                    ignore.case = TRUE)]
OpenTargets_KidneyCancer2Gene_lit_lib <- OpenTargets_KidneyCancer2Gene_lit_lib[grep("adrenal", 
                                                                                    names(OpenTargets_KidneyCancer2Gene_lit_lib), 
                                                                                    ignore.case = TRUE,
                                                                                    invert = TRUE)]
OpenTargets_KidneyCancer2Gene_lit_lib <- OpenTargets_KidneyCancer2Gene_lit_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                                    names(OpenTargets_KidneyCancer2Gene_lit_lib), 
                                                                                    ignore.case = TRUE, 
                                                                                    invert = TRUE)] 
names(OpenTargets_KidneyCancer2Gene_lit_lib) <- paste0(names(OpenTargets_KidneyCancer2Gene_lit_lib), "[OpenTargets_literature]")



# ## Enrichr
# 
# Enrichr_Disease2Gene_GeoDiseaseSig_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Enrichr_Disease2Gene_GeoDiseaseSig_lib.rds")
# 
# Enrichr_GeoKidneyCancerSignatures_Up <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up[grep("kidney|renal",
#                                                                                        names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Up),
#                                                                                        ignore.case = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Up <- Enrichr_GeoKidneyCancerSignatures_Up[grep("cancer|carcinoma|sarcoma",
#                                                                                   names(Enrichr_GeoKidneyCancerSignatures_Up),
#                                                                                   ignore.case = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Up <- Enrichr_GeoKidneyCancerSignatures_Up[grep("adrenal",
#                                                                                   names(Enrichr_GeoKidneyCancerSignatures_Up),
#                                                                                   ignore.case = TRUE,
#                                                                                   invert = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Up <- Enrichr_GeoKidneyCancerSignatures_Up[grep("hereditary|familial|susceptibility|predisposition", 
#                       names(Enrichr_GeoKidneyCancerSignatures_Up), 
#                       ignore.case = TRUE, 
#                       invert = TRUE)] 
# names(Enrichr_GeoKidneyCancerSignatures_Up) <- paste0(names(Enrichr_GeoKidneyCancerSignatures_Up), "[Enrichr_GeoDiseaseSig_Up]")
# 
# Enrichr_GeoKidneyCancerSignatures_Down <- Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down[grep("kidney|renal", names(Enrichr_Disease2Gene_GeoDiseaseSig_lib$Down),
#                                                                                            ignore.case = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Down <- Enrichr_GeoKidneyCancerSignatures_Down[grep("cancer|carcinoma|sarcoma",
#                                                                                       names(Enrichr_GeoKidneyCancerSignatures_Down),
#                                                                                       ignore.case = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Down <- Enrichr_GeoKidneyCancerSignatures_Down[grep("adrenal",
#                                                                                       names(Enrichr_GeoKidneyCancerSignatures_Down),
#                                                                                       ignore.case = TRUE,
#                                                                                       invert = TRUE)]
# Enrichr_GeoKidneyCancerSignatures_Down <- Enrichr_GeoKidneyCancerSignatures_Down[grep("hereditary|familial|susceptibility|predisposition", 
#                       names(Enrichr_GeoKidneyCancerSignatures_Down), 
#                       ignore.case = TRUE, 
#                       invert = TRUE)] 
# names(Enrichr_GeoKidneyCancerSignatures_Down) <- paste0(names(Enrichr_GeoKidneyCancerSignatures_Down), "[Enrichr_GeoDiseaseSig_Down]")



## Intogen

Intogen_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/Intogen_Disease2Gene_lib.rds")
Intogen_KidneyCancer2Gene_lib <- Intogen_Disease2Gene_lib[grep("kidney|renal", names(Intogen_Disease2Gene_lib), 
                                                               ignore.case = TRUE)]
Intogen_KidneyCancer2Gene_lib <- Intogen_KidneyCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                    names(Intogen_KidneyCancer2Gene_lib), 
                                                                    ignore.case = TRUE)]
Intogen_KidneyCancer2Gene_lib <- Intogen_KidneyCancer2Gene_lib[grep("adrenal", 
                                                                    names(Intogen_KidneyCancer2Gene_lib), 
                                                                    ignore.case = TRUE,
                                                                    invert = TRUE)]
Intogen_KidneyCancer2Gene_lib <- Intogen_KidneyCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                    names(Intogen_KidneyCancer2Gene_lib), 
                                                                    ignore.case = TRUE, 
                                                                    invert = TRUE)] 
names(Intogen_KidneyCancer2Gene_lib) <- paste0(names(Intogen_KidneyCancer2Gene_lib), "[Intogen]")



## PharmGKB

PharmGKB_Disease2Gene_lib <- readRDS("InputFiles/Enrichment_analysis_libraries/PharmGKB_Disease2Gene_lib.rds")
PharmGKB_KidneyCancer2Gene_lib <- PharmGKB_Disease2Gene_lib[grep("kidney|renal", 
                                                                 names(PharmGKB_Disease2Gene_lib), 
                                                                 ignore.case = TRUE)]
PharmGKB_KidneyCancer2Gene_lib <- PharmGKB_KidneyCancer2Gene_lib[grep("cancer|carcinoma|sarcoma", 
                                                                      names(PharmGKB_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE)]
PharmGKB_KidneyCancer2Gene_lib <- PharmGKB_KidneyCancer2Gene_lib[grep("adrenal", 
                                                                      names(PharmGKB_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE,
                                                                      invert = TRUE)]
PharmGKB_KidneyCancer2Gene_lib <- PharmGKB_KidneyCancer2Gene_lib[grep("hereditary|familial|susceptibility|predisposition", 
                                                                      names(PharmGKB_KidneyCancer2Gene_lib), 
                                                                      ignore.case = TRUE, 
                                                                      invert = TRUE)] 
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
                                      OpenTargets_KidneyCancer2Gene_RNA_lib, OpenTargets_KidneyCancer2Gene_lit_lib, 
                                      Intogen_KidneyCancer2Gene_lib,
                                      PharmGKB_KidneyCancer2Gene_lib, ThETA_KidneyCancer2Gene_lib)





# Save as RDS file
if(!dir.exists("InputFiles/Enrichment_analysis_libraries/")){dir.create("InputFiles/Enrichment_analysis_libraries/", recursive = TRUE)}
saveRDS(Enrichment_KidneyCancer2Gene_lib, "InputFiles/Enrichment_analysis_libraries/Disease2Gene_KidneyCancer_lib.rds")



print(warnings())