set.seed(5081)



# Script to build protein-protein interaction network from STRING database



# Load libraries
library(biomaRt)
library(igraph)



# Download the protein-protein interactions from stringdb
if(!dir.exists("Databases/StringDB/")){dir.create("Databases/StringDB", recursive = TRUE)}

if(!file.exists("Databases/StringDB/9606.protein.links.detailed.v12.0.txt.gz")){
  download.file(url = "https://stringdb-static.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz",
                destfile = "Databases/StringDB/9606.protein.links.detailed.v12.0.txt.gz", method = "wget")
}


# Read the STRING protein-protein interactions
String_ppi <- read.table("Databases/StringDB/9606.protein.links.detailed.v12.0.txt.gz", header = TRUE)


# Filter all interactions with scores greater than 0 and sourced from database/experiments
String_ppi <- String_ppi[(String_ppi$combined_score > 400), ]


# Get list of proteins in the PPI
String_proteins <- sort(unique(c(String_ppi$protein1, String_ppi$protein2)))


# Create mapping from Ensembl Peptide ID ID to Ensembl Gene ID
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
ensemblPeptideID_2_ensemblGeneID <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id"), 
                                          mart = ensembl, 
                                          filters = "ensembl_peptide_id", 
                                          values = gsub("^9606.", "", String_proteins))
ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id <- paste0("9606.", ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)

String_ppi$Node1_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein1, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]
String_ppi$Node2_ensembl_gene_id <- ensemblPeptideID_2_ensemblGeneID$ensembl_gene_id[match(String_ppi$protein2, ensemblPeptideID_2_ensemblGeneID$ensembl_peptide_id)]


# Extract the STRING PPI
String_ppi_Net <- na.exclude(String_ppi[, c("protein1", "protein2", "Node1_ensembl_gene_id", "Node2_ensembl_gene_id", "combined_score")])
String_ppi_Net$Node1_type <- "gene"
String_ppi_Net$Node2_type <- "gene"
String_ppi_Net$Edge_type <- "undirected"
String_ppi_Net <- String_ppi_Net[, c("Node1_ensembl_gene_id", "Node1_type", "Node2_ensembl_gene_id", "Node2_type", "Edge_type")]

# Convert to graph object
String_ppi_Net <- graph_from_data_frame(String_ppi_Net[, c("Node1_ensembl_gene_id", "Node2_ensembl_gene_id")], directed = FALSE)
String_ppi_Net <- simplify(String_ppi_Net, 
                           remove.loops = TRUE, 
                           remove.multiple	= TRUE) # remove loops and multi-edges


# Find the connected components
ppi_comps <- components(String_ppi_Net)

# Identify the largest component
largest_ppi_comp <- which.max(ppi_comps$csize)

# Extract the largest component
String_ppi_Net <- induced.subgraph(String_ppi_Net, which(ppi_comps$membership == largest_ppi_comp))



print(paste("Network size (vertices, edges):", vcount(String_ppi_Net), ecount(String_ppi_Net)))

if(!dir.exists("InputFiles/Networks/")){dir.create("InputFiles/Networks/", recursive = TRUE)}
saveRDS(String_ppi_Net, "InputFiles/Networks/STRING_PPI_Net.rds")


network_parameters <- data.frame(t(data.frame(nodeNumber = vcount(String_ppi_Net),
                                              edgeNumber = ecount(String_ppi_Net),
                                              diameter = diameter(graph = String_ppi_Net, directed = FALSE),
                                              radius = radius(graph = String_ppi_Net),
                                              averagePathLength = mean_distance(graph = String_ppi_Net),
                                              clusteringCoefficient = transitivity(graph = String_ppi_Net, type = "global", isolates = "zero"),
                                              density = edge_density(graph = String_ppi_Net),
                                              connectedComponents = count_components(graph = String_ppi_Net),
                                              largest_component_size = components(String_ppi_Net)$csize[1],
                                              powerFit = fit_power_law(x = degree(String_ppi_Net)+1)$alpha,
                                              degree = centr_degree(graph = String_ppi_Net)$centralization,
                                              closeness = centr_clo(graph = String_ppi_Net)$centralization,
                                              betweenness = centr_betw(graph = String_ppi_Net)$centralization)))

network_parameters <- tibble::rownames_to_column(network_parameters, "Parameters")
colnames(network_parameters) <- c("Parameters", "Value")
network_parameters$Value <- round(network_parameters$Value, 3)

if(!dir.exists("OutputFiles/Tables/")){
  dir.create("OutputFiles/Tables/", recursive = TRUE)
}
write.csv(network_parameters, "OutputFiles/Tables/STRING_PPI_Net_params.csv", row.names = FALSE)


print(warnings())