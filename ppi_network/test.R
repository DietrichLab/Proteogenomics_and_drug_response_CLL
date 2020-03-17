# tests

library(tidyverse)
library(igraph)

outputFolder <- "output"

# read the graph gml
graph <- read_graph(file = file.path(outputFolder, "CLL_correlation_graph_cutoff_0.5.gml"), format = "gml")

# get the node names
nodes <- V(graph) %>% names()

# check, if they are identical to high variance detection
stopifnot(all(highVarProteins %in% nodes) && identical(length(nodes), length(highVarProteins)))

# get trisomy12 medians for 30 proteins
pois <- sample(nodes, 30)

# get the quant data
df <- wideFormat(mae[ , , "proteomics"], colDataCols = "trisomy12") %>%
  as.data.frame() %>%
  dplyr::filter(trisomy12 == 1) %>%
  dplyr::select(-trisomy12) %>%
  column_to_rownames("primary")
colnames(df) <- gsub("proteomics_", "", colnames(df))

medians.data <- df[, pois] %>%
  as.matrix() %>%
  colMedians(na.rm = TRUE)

# chek for identity
medians.graph <- lapply(X = pois,
                        FUN = function (p) {
                          V(graph)[[p]]$trisomy12
                        }) %>% unlist()

stopifnot(identical(medians.data, medians.graph))
