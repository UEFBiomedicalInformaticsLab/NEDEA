# Create PPI network from CPTAC raw counts

# Load libraries
library(biomaRt)
library(tidyverse)
library(edgeR)
library(DESeq2)
library(WGCNA)
library(igraph)

disease <- "LungCancer"


# Read the raw counts
expr_raw <- read.csv(paste0("Databases/CPTAC/CPTAC_", disease, "_RNAseq_counts.csv"), check.names = FALSE)

# Read clinical data for the samples
sample_annot <- read.csv(paste0("Databases/CPTAC/CPTAC_", disease, "_clinical_data.csv"))
# sample_annot$submitter_id[!sample_annot$submitter_id %in% colnames(expr_raw)]
sample_annot <- column_to_rownames(sample_annot, "submitter_id")
sample_annot <- sample_annot[order(row.names(sample_annot)), ]


# colnames(expr_raw)[!colnames(expr_raw) %in% sample_annot$submitter_id]



# Extract ensemble ids without version
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 102)
ensemblIDver_2_ensemblID <- getBM(attributes = c("ensembl_gene_id",  "ensembl_gene_id_version"), mart = ensembl, filters = "ensembl_gene_id_version", values = expr_raw$gene_id)



# Clean expression matrix
expr_raw$ensemblID <- ensemblIDver_2_ensemblID$ensembl_gene_id[match(expr_raw$gene_id, ensemblIDver_2_ensemblID$ensembl_gene_id_version)]
expr_raw <- expr_raw[!is.na(expr_raw$ensemblID), ] # Remove genes with no matching ensembl ID [mostly removed PAR genes]
expr_raw <- expr_raw[expr_raw$gene_type == "protein_coding", ]
row.names(expr_raw) <- NULL
expr_raw <- column_to_rownames(expr_raw, "ensemblID")
expr_raw <- expr_raw[, !colnames(expr_raw) %in% c("gene_id", "gene_name","gene_type")]
expr_raw <- expr_raw[, order(names(expr_raw))]


# Low count genes filter 
keep <- rowSums(cpm(expr_raw) > 1) > round(0.05 * ncol(expr_raw)) # Filter genes not having 1 CPM in at least 5% of the sample 
expr_raw_filt <- expr_raw[keep,]
dim(expr_raw_filt) 

colnames(expr_raw_filt)[!colnames(expr_raw_filt) %in% row.names(sample_annot)]
row.names(sample_annot)[!row.names(sample_annot) %in% colnames(expr_raw_filt)]

# Normalization
dds <- DESeqDataSetFromMatrix(countData = expr_raw_filt, colData = sample_annot, design = ~ 1) # Check design
# dds <- DESeq(dds)
vsd <- varianceStabilizingTransformation(dds, blind = FALSE) # Check blind
expr_vst <- assay(vsd)



hard_threshold <- pickHardThreshold(t(expr_vst), RsquaredCut = 0.8)

print(paste0("Hard threshold used = ", hard_threshold$cutEstimate))

adjacency_unweight <- signumAdjacencyFunction(cor(t(expr_vst)), hard_threshold$cutEstimate)

network <- graph_from_adjacency_matrix(adjacency_unweight, mode = "undirected")


if(!dir.exists("InputFiles/Networks")){dir.create("InputFiles/Networks", recursive = TRUE)} 
saveRDS(network, paste0("InputFiles/Networks/CPTAC_rawCount_", disease, "_unweighted.rds"))
