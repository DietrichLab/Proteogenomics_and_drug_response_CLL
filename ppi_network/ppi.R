# generate PPI network

library(tidyverse)
library(doMC)
library(mixtools)
library(MultiAssayExperiment)
library(matrixStats)
library(Hmisc)
library(igraph)
library(ggthemes)

sourceDir <- function(path, ...) {
  for (nm in list.files(path, pattern = "\\.[RrSsQq]$")) {
    source(file.path(path, nm), ...)
  }
}
sourceDir(path = "functions")

outputFolder <- "output"


# load proteomics data
load("multiomics_MAE.RData")
mae <- multiomics_MAE
rm(multiomics_MAE)

# get the clustering
# get clustering
load("CLL_Proteomics_ConsensusClustering.RData")
colData(mae)$PG <- as.factor(CCP_group5[rownames(colData(mae))])
rm(prot_cor_meta, CCP_group5, CCP_group6_RNA)


# extract protein expressions
proteins.wide <- mae[, , "proteomics"] %>%
  wideFormat() %>%
  as.data.frame() %>%
  column_to_rownames("primary")
colnames(proteins.wide) <- gsub("proteomics_", "", colnames(proteins.wide))


# log2 expression distribution QC
proteins.wide %>%
  t() %>%
  as.data.frame() %>%
  gather(sample, value) %>%
  ggplot(aes(x = value, group = sample, color = sample)) +
  geom_density() +
  xlim(c(-2.5, 2.5)) +
  xlab("log2 expression") +
  theme_tufte() +
  theme(legend.position = "none")


# detect high variance proteins
highVarProteins <- detectHighVarianceProteins(expressions = t(proteins.wide),
                                              nCores = 2,
                                              maxit = 1000,
                                              outputFolder = file.path(outputFolder, "high_variance_detection"))

proteins.highVar <- proteins.wide[, highVarProteins]

# export labeled proteins for visualization
proteins.wide %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column("gene_symbol") %>%
  mutate(highvar = gene_symbol %in% highVarProteins) %>%
  select(gene_symbol, highvar, everything()) %>%
  write.table(file = "proteins.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# correlate expressions
cor.obj <- rcorr(as.matrix(proteins.highVar), type = "pearson")
cor.mat <- cor.obj$r

# cluster the correlation matrix
clust <- hclust(d = as.dist(1 - cor.obj$r), method = "ward.D2")
cor.mat.clust <- cor.obj$r[clust$order, clust$order]
dim(cor.mat.clust)

# look at r distribution
cor.mat.qq <- cor.mat[upper.tri(cor.mat)]

pdf(file = file.path(outputFolder, "correlation_cutoff_qq.pdf"))
qqnorm(cor.mat.qq)
qqline(cor.mat.qq)
# 0.5 seems to be a good cutoff for the network generation
cutoff <- 0.5
abline(h = cutoff)
dev.off()


# build a network
cor.mat.g <- cor.mat
cor.mat.g[cor.mat.g < cutoff] <- 0
g <- graph_from_adjacency_matrix(adjmatrix = cor.mat.g,
                                 mode = "upper",
                                 weighted = TRUE,
                                 diag = FALSE,
                                 add.colnames = NULL)

stopifnot(identical(colnames(cor.mat), names(V(g))))
stopifnot(identical(colnames(cor.mat), highVarProteins))
stopifnot(identical(names(V(g)), highVarProteins))


# get some expression medians
g <- storeGroupMediansInGraph(mae = mae,
                              proteins = highVarProteins,
                              g = g,
                              type = "trisomy12",
                              required = 1)
g <- storeGroupMediansInGraph(mae = mae,
                              proteins = highVarProteins,
                              g = g,
                              type = "IGHV",
                              required = "M")
g <- storeGroupMediansInGraph(mae = mae,
                              proteins = highVarProteins,
                              g = g,
                              type = "TP53",
                              required = 1)

for (i in 1:6) {
  g <- storeGroupMediansInGraph(mae = mae,
                                proteins = highVarProteins,
                                g = g,
                                type = "PG",
                                required = i)
}


# save the graph as gml
write_graph(graph = g,
            file = file.path(outputFolder, paste0("CLL_correlation_graph_cutoff_", cutoff, ".gml")),
            format = "gml")


# write correlation matrix to file
write.table(x = cor.mat.clust,
            file = file.path(outputFolder, "CLL_correlation_matrix_full.txt"),
            quote = FALSE,
            sep = "\t",
            row.names = TRUE,
            col.names = NA)
