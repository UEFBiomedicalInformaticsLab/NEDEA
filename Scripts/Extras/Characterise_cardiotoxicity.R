library(igraph)


# Read the drug type information
DrugBank_drug_type <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_drug_type <- DrugBank_drug_type$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


# Read drug target interactions
drug_target_ixn <- readRDS("InputFiles/Associations/DrugBank_Drug_Target_associations.rds")



# Read the DDI data
DrugBank_ddi <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")
DrugBank_ddi <- DrugBank_ddi$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,4)] <- c("Drug1_DrugBank_id", "Drug2_DrugBank_id")
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% DrugBank_drug_type$primary_key & DrugBank_ddi$Drug2_DrugBank_id %in% DrugBank_drug_type$primary_key, ] # retain only small molecular drugs
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% drug_target_ixn$drugbank_drug_id & DrugBank_ddi$Drug2_DrugBank_id %in% drug_target_ixn$drugbank_drug_id, ] # retain drugs with known targets



# Extract DDI with toxicity
DrugBank_ddi <- DrugBank_ddi[grep("cardiotoxicity|cardiotoxic activities", DrugBank_ddi$description, ignore.case = TRUE), ]

tmp1 <- gsub(pattern = ".*(risk or severity of .+toxicity can be increased).*|.*(risk or severity of .+toxicity can be decreased).*|.*(increase the .+toxic activities).*|.*(decrease the .+toxic activities).*",
             replacement = "\\1\\2\\3\\4",
             x = DrugBank_ddi$description,
             ignore.case = TRUE)

sort(table(tmp1), decreasing = TRUE)





# Separate the DDI that increase and decrease the toxicity
terms_increase <- unique(tmp1)[grep("increase", unique(tmp1), ignore.case = TRUE)]
terms_decrease <- unique(tmp1)[grep("decrease", unique(tmp1), ignore.case = TRUE)]

DrugBank_increase_ddi <- DrugBank_ddi[grep(pattern = paste(terms_increase, collapse = "|"), DrugBank_ddi$description),]
DrugBank_decrease_ddi <- DrugBank_ddi[grep(pattern = paste(terms_decrease, collapse = "|"), DrugBank_ddi$description),]

rm(tmp1)







# Create a network using the possible target pairs


DDI_increase_net <- data.frame()

for(i in 1:nrow(DrugBank_increase_ddi)){
  drug1 <- DrugBank_increase_ddi[i,"Drug1_DrugBank_id", drop = TRUE]
  drug2 <- DrugBank_increase_ddi[i,"Drug2_DrugBank_id", drop = TRUE]
  
  if(drug1 == drug2){ print(i) }
  drug1_targets <- unique(drug_target_ixn[drug_target_ixn$drugbank_drug_id %in% drug1, "ensembl_gene_id", drop = TRUE])
  drug2_targets <- unique(drug_target_ixn[drug_target_ixn$drugbank_drug_id %in% drug2, "ensembl_gene_id", drop = TRUE])
  
  
  # Generate unique combinations
  tmp1 <- expand.grid(drug2_targets, drug2_targets)
  colnames(tmp1) <- c("target1", "target2")
  tmp1 <- tmp1[!duplicated(t(apply(tmp1, 1, sort))), ]
  
  tmp1$comb_name <- paste(drug1, drug2, sep = "_")
  DDI_increase_net <- rbind(DDI_increase_net, tmp1)
}

DDI_increase_net <- DDI_increase_net %>% count(target1, target2, sort = TRUE)
DDI_increase_net <- graph.data.frame(DDI_increase_net, directed = FALSE)
createNetworkFromIgraph(igraph = DDI_increase_net)




DDI_decrease_net <- data.frame()

for(i in 1:nrow(DrugBank_decrease_ddi)){
  drug1 <- DrugBank_decrease_ddi[i,"Drug1_DrugBank_id", drop = TRUE]
  drug2 <- DrugBank_decrease_ddi[i,"Drug2_DrugBank_id", drop = TRUE]
  
  if(drug1 == drug2){ print(i) }
  drug1_targets <- unique(drug_target_ixn[drug_target_ixn$drugbank_drug_id %in% drug1, "ensembl_gene_id", drop = TRUE])
  drug2_targets <- unique(drug_target_ixn[drug_target_ixn$drugbank_drug_id %in% drug2, "ensembl_gene_id", drop = TRUE])
  

  # Generate unique combinations
  tmp1 <- expand.grid(drug2_targets, drug2_targets)
  colnames(tmp1) <- c("target1", "target2")
  tmp1 <- tmp1[!duplicated(t(apply(tmp1, 1, sort))), ]
  
  tmp1$comb_name <- paste(drug1, drug2, sep = "_")
  DDI_decrease_net <- rbind(DDI_decrease_net, tmp1)
}

DDI_decrease_net <- DDI_decrease_net %>% count(target1, target2, sort = TRUE)
DDI_decrease_net <- graph.data.frame(DDI_decrease_net, directed = FALSE)
createNetworkFromIgraph(igraph = DDI_decrease_net)

