library(biomaRt)


ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
peptide_IDs <- getBM(attributes = c("ensembl_gene_id", "ensembl_peptide_id", "gene_biotype"), 
                                          mart = ensembl)
peptide_IDs <- peptide_IDs[peptide_IDs$gene_biotype %in% "protein_coding", ]
peptide_IDs <- peptide_IDs[peptide_IDs$ensembl_peptide_id != "", ]



peptide_sequences <- getSequence(id = peptide_IDs$ensembl_peptide_id, 
                                 type = "ensembl_peptide_id", 
                                 seqType = "peptide", 
                                 mart = ensembl)