set.seed(5081)



# Script to build protein-protein interaction network from STRING based on information scores of databases



# Load libraries
library(igraph)
library(biomaRt)





if(!dir.exists("Databases/StringDB/")){dir.create("Databases/StringDB", recursive = TRUE)}

if(!file.exists("Databases/StringDB/9606.protein.links.detailed.v11.5.txt.gz")){
  download.file(url = "https://stringdb-static.org/download/protein.links.detailed.v11.5/9606.protein.links.detailed.v11.5.txt.gz",
                destfile = "Databases/StringDB/9606.protein.links.detailed.v11.5.txt.gz", method = "wget")
}

# if(!file.exists("Databases/StringDB/9606.protein.aliases.v11.5.txt.gz")){
#   download.file(url = "https://stringdb-static.org/download/protein.aliases.v11.5/9606.protein.aliases.v11.5.txt.gz",
#                 destfile = "Databases/StringDB/9606.protein.aliases.v11.5.txt.gz", method = "wget")
# }



# # Read protein Id mappings
# String_proteins <- read.table("Databases/StringDB/9606.protein.aliases.v11.5.txt.gz", fill = TRUE, sep = "\t")
# colnames(String_proteins) <- c("string_protein_id", "alias", "source")
# String_proteins <- String_proteins[String_proteins$source == "Ensembl_HGNC_Ensembl_ID(supplied_by_Ensembl)",]



# Read all PPI in String DB along with scores
String_ppi <- read.table("Databases/StringDB/9606.protein.links.detailed.v11.5.txt.gz", header = TRUE)


# Filter all interactions with scores greater than 0 and sourced from database/experimets
String_ppi <- String_ppi[(String_ppi$database > 0 | String_ppi$experimental > 0), ] #

String_proteins <- sort(unique(c(String_ppi$protein1, String_ppi$protein2)))

# Create mapping from Ensembl Peptide ID ID to Ensembl Gene ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensemblPeptideID_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                                          mart = ensembl, 
                                          filters = "ensembl_peptide_id", 
                                          values = gsub("^9606.", "", String_proteins))
ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id <- paste0("9606.", ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)


# String_ppi$Node1_ensembl_gene_id <- String_proteins$alias[match(String_ppi$protein1, String_proteins$string_protein_id)]
# String_ppi$Node2_ensembl_gene_id <- String_proteins$alias[match(String_ppi$protein2, String_proteins$string_protein_id)]

String_ppi$Node1_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein1, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]
String_ppi$Node2_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein2, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]


String_ppi_Net <- na.exclude(String_ppi[, c("protein1", "protein2", "Node1_ensembl_gene_id", "Node2_ensembl_gene_id")])
String_ppi_Net$Node1_type <- "gene"
String_ppi_Net$Node2_type <- "gene"
String_ppi_Net$Edge_type <- "undirected"
String_ppi_Net <- String_ppi_Net[, c("Node1_ensembl_gene_id", "Node1_type", "Node2_ensembl_gene_id", "Node2_type", "Edge_type")]

# Convert to graph object
String_ppi_Net <- graph_from_data_frame(String_ppi_Net[, c("Node1_ensembl_gene_id", "Node2_ensembl_gene_id")], directed = FALSE)
String_ppi_Net <- simplify(String_ppi_Net, remove.loops = TRUE, remove.multiple	= TRUE) # remove loops and multi-edges
print(paste("Network size (vertices, edges):", vcount(String_ppi_Net), ecount(String_ppi_Net)))

if(!dir.exists("InputFiles/Networks/")){dir.create("InputFiles/Networks/", recursive = TRUE)}
saveRDS(String_ppi_Net, "InputFiles/Networks/STRING_PPI_Net_database.rds")



print(warnings())