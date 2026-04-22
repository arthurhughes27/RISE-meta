aggregate_to_geneset = function(df, genesets, geneset_names, FUN = mean) {
  df <- as.data.frame(df)
  
  # ---- no aggregation case ----
  if (is.character(FUN) && FUN == "none") {
    return(df)
  }
  
  geneset_names <- as.character(geneset_names)
  
  out <- imap_dfc(genesets, function(genes, i) {
    present <- intersect(genes, colnames(df))
    if (length(present) == 0) {
      return(NULL)
    }
    
    mat <- df[, present, drop = FALSE]
    vec <- apply(mat, 1, function(r)
      FUN(r, na.rm = TRUE))
    
    # Replace NaN (all-NA rows) with proper NA
    vec[is.nan(vec)] <- NA_real_
    
    tibble::tibble(!!geneset_names[[i]] := vec)
  })
  
  as.data.frame(out, check.names = FALSE)
}