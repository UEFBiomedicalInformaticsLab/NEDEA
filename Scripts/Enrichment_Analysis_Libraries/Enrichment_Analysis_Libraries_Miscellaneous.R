set.seed(5081)
rm(list = ls())



# Libraries for enrichment analysis of miscellaneous gene sets





# Load libraries
library(unixtools)
library(org.Hs.eg.db)
library(openxlsx)
# library(msigdbr)
# library(httr)



# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")



# Create gene ID mappings
entrezId_2_ensemblId <- as.data.frame(org.Hs.egENSEMBL)
entrezId_2_geneSymbol <- as.data.frame(org.Hs.egSYMBOL)
geneSymbol_2_ensemblId <- merge(entrezId_2_ensemblId, entrezId_2_geneSymbol, by = "gene_id")





## Enrichment of drugable genome list from DGIdb -------------------------------------------------------------

# Retrieve the gene categories
if(!dir.exists("Databases/DGIdb")){dir.create("Databases/DGIdb", recursive = TRUE)}
if(!file.exists("Databases/DGIdb/categories.tsv")){
  download.file(url = "https://www.dgidb.org/data/monthly_tsvs/2022-Feb/categories.tsv",
                destfile = "Databases/DGIdb/categories.tsv", method = "wget")
}


dgidb_drugable_genome <- read.table("Databases/DGIdb/categories.tsv", sep = "\t", header = TRUE, quote = "")
dgidb_drugable_genome <- dgidb_drugable_genome[dgidb_drugable_genome$category == "DRUGGABLE GENOME", ]
dgidb_drugable_genome$ensembl_gene_id <- geneSymbol_2_ensemblId$ensembl_id[match(dgidb_drugable_genome$entrez_gene_symbol, geneSymbol_2_ensemblId$symbol)]
dgidb_drugable_genome <- dgidb_drugable_genome$ensembl_gene_id[!is.na(dgidb_drugable_genome$ensembl_gene_id)]





## Enrichment of Very Important Pharmacogenes (VIP) list from PharmGKB -------------------------------------------------------------

if(!dir.exists("Databases/PharmGKB/")){dir.create("Databases/PharmGKB/", recursive = TRUE)}
if(!file.exists("Databases/PharmGKB/genes.tsv")){
  download.file(url = "https://api.pharmgkb.org/v1/download/file/data/genes.zip",
                destfile = "Databases/PharmGKB/genes.zip", method = "wget")
  unzip("Databases/PharmGKB/genes.zip", exdir = "Databases/PharmGKB/", file = "genes.tsv")
}


PharmGKB_genes <- read.table("Databases/PharmGKB/genes.tsv", sep = "\t", fill = TRUE, quote = "", header = TRUE)
PharmGKB_vip <- PharmGKB_genes[PharmGKB_genes$Is.VIP == "Yes", "Ensembl.Id"]
PharmGKB_vip <- PharmGKB_vip[(PharmGKB_vip != "")]
PharmGKB_vip <- PharmGKB_vip[grep("\",", PharmGKB_vip, invert = TRUE)]





## Enrichment of genes related to ageing from GenAge db -------------------------------------------------------------

if(!dir.exists("Databases/GenAge/")){dir.create("Databases/GenAge/", recursive = TRUE)}
if(!file.exists("Databases/GenAge/genage_human.csv")){
  download.file(url = "https://genomics.senescence.info/genes/human_genes.zip",
                destfile = "Databases/GenAge/human_genes.zip", method = "wget")
  unzip("Databases/GenAge/human_genes.zip", exdir = "Databases/GenAge/", file = "genage_human.csv")
}

GenAge_genes <- read.csv("Databases/GenAge/genage_human.csv")
GenAge_genes$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(GenAge_genes$entrez.gene.id, entrezId_2_ensemblId$gene_id)]
GenAge_genes <- as.vector(na.exclude(GenAge_genes$ensembl_gene_id))





## Enrichment of hallmarkers of cancer genes from HOC db -------------------------------------------------------------
# Ref: https://doi.org/10.1093/database/baaa045
if(!dir.exists("Databases/HOCdb/")){dir.create("Databases/HOCdb/", recursive = TRUE)}

# Manually download supplementary file from https://doi.org/10.1093/database/baaa045
if(!file.exists("Databases/HOCdb/suppl_data_baaa045.zip")){
    warning(paste0("Download supplementary file: Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and Pathway.xlsx"))
}
unzip("Databases/HOCdb/suppl_data_baaa045.zip", exdir = "Databases/HOCdb/", file = "Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and Pathway.xlsx")

HOC_genes <- read.xlsx("Databases/HOCdb/Halifax-curation.Table S2. Hallmark of Cancer Data - Gene and Pathway.xlsx", sheet = "HP-gene")
HOC_genes$ensembl_gene_id <- entrezId_2_ensemblId$ensembl_id[match(HOC_genes$ID, entrezId_2_ensemblId$gene_id)]
HOC_genes <- unique(HOC_genes$ensembl_gene_id)





# Compile all miscellaneous gene list into one list
miscellaneous_gene_lib <- list(DGdb_DrugableGenome = dgidb_drugable_genome,
                               PharmGKB_vip = PharmGKB_vip,
                               GenAge_genes = GenAge_genes,
                               Hallmaker_of_Cancer_genes = HOC_genes)

if(!dir.exists("InputFiles/Enrichment_Analysis_Libraries/")){dir.create("InputFiles/Enrichment_Analysis_Libraries/", recursive = TRUE)}
saveRDS(miscellaneous_gene_lib, "InputFiles/Enrichment_Analysis_Libraries/miscellaneous_gene_lib.rds")


print(warnings())