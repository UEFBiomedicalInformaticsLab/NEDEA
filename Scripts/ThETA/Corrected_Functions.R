# Functions corrected from TheTA


tissue.specific.scores <- function(disease_genes, ppi_network, directed_network = F,
                                   tissue_expr_data,  dis_relevant_tissues, W, cutoff = 1.6,
                                   selected_tissues = NULL, verbose = FALSE) {
  if(is.null(disease_genes) | is.null(tissue_expr_data) | is.null(ppi_network) | 
     (is.null(dis_relevant_tissues) & is.null(selected_tissues)))
    stop('Incomplete data input.')
  # compile the wsp
  if(is.null(selected_tissues) & !is.null(dis_relevant_tissues)) {
    wsp_list = weighted.shortest.path(disease_genes, ppi_network,
                                      directed_network, tissue_expr_data,
                                      dis_relevant_tissues, W, cutoff, selected_tissues = NULL, verbose)    # Changed tissues to selected_tissues
  }
  else if(is.character(selected_tissues)){
    print("New function..")
    selected_tissues <- intersect(colnames(tissue_expr_data), selected_tissues)
    if(length(selected_tissues) == 0) stop('Please provide a set of valid tissue names.')
    if(verbose)
      print(paste("Number of selected tissues:", length(selected_tissues), sep=""))
    wsp_list = weighted.shortest.path(disease_genes = disease_genes, 
                                      ppi_network = ppi_network,
                                      directed_network = directed_network, 
                                      tissue_expr_data = tissue_expr_data,
                                      dis_relevant_tissues = dis_relevant_tissues, 
                                      W = W, cutoff = NULL, 
                                      selected_tissues = selected_tissues, verbose)
  }
  else 
    print("Please, provide a valid list of tissue names.")
  
  tissue_scores <- apply(wsp_list$shortest_paths, 2, function(x) {
    q <- stats::quantile(x)
    outliers <- q[4]+(1.5*(q[4]-q[2]))
    ts <- scales::rescale(x, from=c(min(x), outliers), to=c(1,0))
    ts[ts< 0] <- 0
    ts
  })
  dt <- data.frame(tissue_scores, avg_tissue_score = rowMeans(tissue_scores[,,drop=F]))
  return(dt)
}