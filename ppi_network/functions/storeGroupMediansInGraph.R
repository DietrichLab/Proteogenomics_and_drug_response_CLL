storeGroupMediansInGraph <- function (mae, proteins, g, type, required = 1) {
  
  df <- wideFormat(mae[ , , "proteomics"], colDataCols = type) %>%
    as.data.frame() %>%
    dplyr::filter(UQ(as.symbol(type)) == required) %>%
    dplyr::select(-UQ(as.symbol(type))) %>%
    column_to_rownames("primary")
  colnames(df) <- gsub("proteomics_", "", colnames(df))
  
  medians <- df[, proteins] %>%
    as.matrix() %>%
    colMedians()
  names(medians) <- colnames(df[, proteins])
  
  g.updated <- set_vertex_attr(graph = g,
                               name = paste0(type, "_", required),
                               index = names(medians),
                               value = medians)
  
  return(g.updated)
}
