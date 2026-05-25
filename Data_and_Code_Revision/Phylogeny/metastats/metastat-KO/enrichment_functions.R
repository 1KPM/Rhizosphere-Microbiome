perform_taxa_enrichment <- function(tax_level, diff_data, tax_data, 
                                    universe = NULL, 
                                    p_cutoff = 0.05, min_counts = 1) {
  merged_data <- diff_data %>% 
    select(FeatureID, Significant, !!sym(tax_level)) %>%
    na.omit()
  
  if(is.null(universe)) {
    universe <- unique(tax_data$FeatureID)
  }
  
  N <- length(universe)
  k <- length(unique(merged_data$FeatureID))
  tax_counts <- table(merged_data[[tax_level]])
  
  valid_taxa <- names(tax_counts)[tax_counts >= min_counts]
  
  hyper_results <- map_df(valid_taxa, function(taxon) {
    M <- sum(tax_data[[tax_level]] == taxon & 
               tax_data$FeatureID %in% universe)
    x <- sum(merged_data[[tax_level]] == taxon & 
               merged_data$Significant)
    
    p_val <- phyper(q = x - 1, m = M, n = N - M, k = k, 
                    lower.tail = FALSE)
    
    gene_ratio <- x / k
    bg_ratio <- M / N
    enrichment_factor <- gene_ratio / bg_ratio
    
    data.frame(
      Taxon = taxon,
      Count = x,
      Total = M,
      GeneRatio = paste(x, k, sep="/"),
      BgRatio = paste(M, N, sep="/"),
      EnrichmentFactor = enrichment_factor,
      pvalue = p_val
    )
  })
  
  hyper_results$FDR <- p.adjust(hyper_results$pvalue, method = "fdr")
  
  hyper_results %>%
    mutate(
      TaxonomicLevel = tax_level,
      Significant = FDR < p_cutoff
    )
}