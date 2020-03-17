# CLL_acetylation


# Libraries
library(missMDA)
library(Rtsne)
library(data.table)
library(AnnotationHub)
library(ensembldb)
library(ggplot2)
library(ggrepel)
library(dplyr)
library(limma)
library(gridExtra)
library(grid)
library(cowplot)


## Functions
abs_axis <- function(i) {
  abs(i)
}


## Data
# Ensembl
hub = AnnotationHub()
# dm <- query(hub, c("Ensembl", "Homo sapiens"))
# dm <- data.frame(id = dm$ah_id, title = dm$title)
edb <- hub[["AH60977"]]


# SetIDs
col_list <- c("#ffbe44", "#00b35b", "#9b2193", "#ec3e6d", "#9997ff", "#753300")
names(col_list) <- paste0('PG', 1:6)

samples_annot <- read.table('~/Downloads/CLL_project/data/patIDmeta_translation_table.txt', header = TRUE)


# Correct group
meta <- read.table('~/Downloads/CLL_project/data/CLL_Meta_Proteomics_2019-12-09.txt', sep = '\t', header = TRUE)

# Change according to the new cluster
samples_annot$Cluster[match(meta$Proteomics_Pat_ID, samples_annot$Labels_proteomics)] = meta$CCP2


## Peptide table
peptab <- fread('~/Downloads/CLL_project/data/Acetylation_search/peptides_table.txt')

# Map set to samples
# Differentially expressed transcripts
# Correct triplicates

tocorrect <- names(which(table(samples_annot$Labels_proteomics) > 1))
samples_annot$is_replicate[samples_annot$Labels_proteomics %in% tocorrect] = TRUE 

tochange <- grep('POOL', colnames(peptab))

set_names <- paste0('Set', rep(c('A', 'B','C', 'D', 'E', 'F', 'G', 'H'), each = 10))
colnames(peptab)[tochange] <- gsub('^.*Set[A-Z]_([^ ])','\\1', colnames(peptab)[tochange])
colnames(peptab)[tochange]  <- paste0(set_names, '_', colnames(peptab)[tochange])


# Quantification values
peptab <- as.data.frame(peptab)
peptabAnnot <- peptab[, c("Peptide sequence", "Protein(s)", "Gene(s)"), ]


gene_ids <- unique(unlist(strsplit(peptabAnnot$`Gene(s)`, ';')))

geneTab <- genes(edb, filter = ~gene_id %in%  gene_ids)
matching_tab <- data.frame('Gene' = geneTab$gene_id, 
                           'Gene_symbol' = geneTab$symbol)


peptabquant <- peptab[, colnames(peptab)[tochange]]
rownames(peptabquant) <- peptab$`Peptide sequence`

total_sample <- ncol(peptabquant)

# Merge replicate
replicates1 <- c('SetA_tmt10plex_130C', 'SetB_tmt10plex_126', 'SetE_tmt10plex_128N')
replicates2 <- c('SetF_tmt10plex_127N', 'SetG_tmt10plex_129C', 'SetH_tmt10plex_131')


replicates1_Tab <- peptabquant[ , replicates1]
merged_rep1 <- apply(replicates1_Tab, 1, mean, na.rm = TRUE)
merged_rep1[is.nan(merged_rep1)] <- NA

replicates2_Tab <- peptabquant[ , replicates2]
merged_rep2 <- apply(replicates2_Tab, 1, mean, na.rm = TRUE)
merged_rep2[is.nan(merged_rep2)] <- NA

# Merge with replicates
non_replicates <- peptabquant[, as.character(samples_annot$TMT.channel[!samples_annot$is_replicate])]

non_replicates <- cbind.data.frame(non_replicates, merged_rep1, merged_rep2)

idChange <- c(ncol(non_replicates)-1, ncol(non_replicates))
colnames(non_replicates)[-idChange] <- 
  as.character(samples_annot$Labels_proteomics[match(colnames(non_replicates)[-idChange], samples_annot$TMT.channel)])

colnames(non_replicates)[idChange] <- c('P0610', 'P0718')



## Acetylation pattern
acetyl_pattern <- '+229.163-187.152'
idx_acetyl <-  grep(pattern = acetyl_pattern, peptab$`Peptide sequence`)


# Reorder
samples_annot <- samples_annot[match(colnames(non_replicates), samples_annot$Labels_proteomics), ]
non_replicates_acetyl <- non_replicates[idx_acetyl, ]


## PCA,TSNE plot
tsne.pca <- function(dat) {
  
  dat <- dat[apply(dat, 1, function(i) sum(is.na(i)) < 0.5 * ncol(dat)), ]
  
  
  ## Impute
  imputedDat <- imputePCA(dat, ncp = 2, scale = TRUE, method = c("EM"), 
                          row.w = NULL, ind.sup=NULL,quanti.sup=NULL,quali.sup=NULL,
                          coeff.ridge = 1, threshold = 1e-06, seed = NULL, nb.init = 1,  
                          maxiter = 1000)
  
  
  # PCA 
  pr = prcomp(t(imputedDat$completeObs), scale. = TRUE)
  
  pr_explained <- round((pr$sdev^2/sum(pr$sdev^2))[c(1,2)], 2) * 100 
  
  
  pca_dat <- data.frame('pca' = pr$x[, c(1,2)], 
                        'samples' = colnames(dat), 
                        'cluster' =  paste0('PG', samples_annot$Cluster[match(colnames(dat), samples_annot$Labels_proteomics)]))
  
  
  
  p <- ggplot(data = pca_dat, aes(x = pca.PC1, y = pca.PC2, col = as.factor(cluster))) + 
    geom_point() + 
    # geom_label_repel(aes(label = samples), size = 2) +
    scale_color_manual(name = '', values = col_list) +
    labs(x = paste0('PC1', ' (',paste0(pr_explained[1], '% of variance explained)')), 
         y = paste0('PC2', ' (',paste0(pr_explained[2], '% of variance explained)')),
         title = 'Acetylation') + 
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5)) + 
    guides(col = guide_legend(override.aes = list(fill = col_list, size = 2)))
  
  
  
  pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_proteomics_PCA_acetylation.pdf'), width = 6, height = 5)
  plot(p)
  dev.off()
  
  
  ## TSNE
  set.seed(123454)
  tsne_out = Rtsne(t(imputedDat$completeObs),
                   # perplexity = round(sqrt(ncol(imputedDat$completeObs)), 2), 
                   perplexity = 30/3,
                   dims = 2)
  
  tsne_dat <- data.frame('tsne' = tsne_out$Y,
                         
                         'samples' = colnames(dat), 
                         'cluster' = paste0('PG', samples_annot$Cluster[match(colnames(dat), samples_annot$Labels_proteomics)]))
  
  p <- ggplot(data = tsne_dat, aes(x = tsne.1, y = tsne.2,  col = as.factor(cluster))) + 
    geom_point() + 
    # geom_label_repel(aes(label = samples), size = 2) + 
    scale_color_manual(name = '', values = col_list) +
    labs(x = 'TSNE-1', y = 'TSNE-2',  title = 'Acetylation') + 
    guides(col = guide_legend(override.aes = list(fill = col_list, size = 2))) +
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5))
  
  
  
  pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_proteomics_TSNE_acetylation.pdf'), width = 6, height = 5)
  plot(p)
  dev.off()
  
  # return(tsne_dat)
}

# tsne.pca(dat = non_replicates_acetyl)
###################################################################################################

### Limma differential abundance 
clusters <-  as.factor(paste0('PG', samples_annot$Cluster))


design_mat <- model.matrix(~0 + clusters)
colnames(design_mat) <- levels(clusters)

# Fit linear model
## Filter by quantification in at least three sets
non_replicates_acetyl <- non_replicates_acetyl[apply(non_replicates_acetyl, 1 , function(i) sum(!is.na(i)) > 18), ]

tocompare <- c("PG1-(PG2+PG3+PG4+PG5+PG6)/5",
               "PG2-(PG1+PG3+PG4+PG5+PG6)/5",
               "PG3-(PG1+PG2+PG4+PG5+PG6)/5",
               "PG4-(PG1+PG2+PG3+PG5+PG6)/5",
               "PG5-(PG1+PG2+PG3+PG4+PG6)/5",
               "PG6-(PG1+PG2+PG3+PG4+PG5)/5")


### Normalize by relative protein abundance at the protein level
# dat_summarized <- lapply(list(non_replicates_acetyl, non_replicates), function(dat) {
#   
#   dat$GeneId <- peptabAnnot$`Gene(s)`[match(rownames(dat), peptabAnnot$`Peptide sequence`)]
#   res_summarized <- as.data.frame(dat %>% group_by(GeneId) %>% summarise_all(median, na.rm = TRUE))
#   rownames(res_summarized) <- res_summarized$GeneId
#   res_summarized$GeneId <- NULL
#   res_summarized
# })
# 
# non_replicates_acetyl_summarized <- dat_summarized[[1]]
# non_replicates_summarized <- dat_summarized[[2]]
# 
# common_genes <- intersect(rownames(non_replicates_acetyl_summarized), rownames(non_replicates_summarized))
# 
# non_replicates_acetyl_summarized <- non_replicates_acetyl_summarized[common_genes, ]
# non_replicates_summarized <- non_replicates_summarized[common_genes, ]
# 
# 
# # Divide acetylated by total relative quantification
# normalized_summarized <- non_replicates_acetyl_summarized - non_replicates_summarized

tofit <- non_replicates_acetyl
fit <- lmFit(tofit, design_mat)

contrasts <- makeContrasts(contrasts = tocompare, levels =  clusters)

fit2 <- contrasts.fit(fit, contrasts)
fit3 <- eBayes(fit2)


# Per gene
resLimma = do.call(rbind, lapply(1:6, function(i) {
  # print(i)
  res <- topTable(fit3, coef = i, number = Inf)
  res <- res[complete.cases(res), ]
  res$Cluster <- as.factor(paste0('PG', i))
  res$Peptide <- rownames(res)
  res$Gene <- peptabAnnot$`Gene(s)`[match(rownames(res), peptabAnnot$`Peptide sequence`)]
  # res$Protein <- peptabAnnot$`Protein(s)`[match(rownames(res), peptabAnnot$`Gene(s)`)]
  res$Protein <- peptabAnnot$`Protein(s)`[match(rownames(res),peptabAnnot$`Peptide sequence`)]
  res$GS <- sapply(res$Gene, function(i) {
    genesToMatch <- unlist(strsplit(i, ';'))
    paste(matching_tab$Gene_symbol[match(genesToMatch, matching_tab$Gene)], collapse =';')
  })
  
  res
}))

sig_level <- c(0.1, 0.05, 0.01)
aggregated_plots <- lapply(sig_level, function(i) {
  
  # print(i)
  countSig <- as.data.frame(resLimma %>% group_by(Cluster, .drop=FALSE) %>% dplyr::filter(adj.P.Val < i) %>% dplyr::count())
  
  countNonSig <- as.data.frame(resLimma %>% group_by(Cluster, .drop=FALSE) %>% dplyr::filter(adj.P.Val > i) %>% dplyr::count())
  
  resp <- sapply(1:6, function(j) {
    
    tab <- matrix(c(countSig$n[j],sum(countSig$n) - countSig$n[j], 
                    countNonSig$n[j],sum(countNonSig$n) - countNonSig$n[j]),
                  byrow = TRUE,
                  nrow = 2,
                  dimnames = list(FC = c("Mod", "NonMod"),
                                  Group = c('This', 'Rest')))
    res <- fisher.test(tab, alternative = 'greater')
    res$p.value
    
  })
  resp <- p.adjust(resp, method = 'BH')
  
  
  
  
  list('CountSig' = countSig, 'FisherP' =  resp)
})


names(aggregated_plots) <- paste0('FDR - ', sig_level)

toplot <- reshape2::melt(lapply(aggregated_plots, function(i) i[[1]]))
fisherTab <-  reshape2::melt(lapply(aggregated_plots, function(i) i[[2]]))
fisherTab$Cluster <- rep(paste0('PG', 1:6), 3)
fisherTab$Sig <- ifelse(fisherTab$value < 0.001, '***', ifelse(fisherTab$value < 0.01, '**', ifelse(fisherTab$value < 0.1, '*', '')))
fisherTab$value <- toplot$value


levToChoose <-  0.1


p1 <- ggplot(subset(toplot, L1 ==  paste0('FDR - ', levToChoose)), aes(x = Cluster, y = value, fill = Cluster)) + 
  geom_bar(stat = 'identity') + 
  geom_text(subset(fisherTab, L1 ==  paste0('FDR - ', levToChoose)),  mapping = aes(label  = Sig), size = 5) +
  scale_fill_manual(values  = col_list) +
  # facet_grid(.~L1) +
  labs(x = '', y = 'Differentially abundant peptides \n(Adjusted p-value < 0.1)', title = 'Acetylation')+ 
  scale_y_continuous(expand = c(0,0, 0.1, 0), breaks = 10 * seq(0, round(max(toplot$value)/10, 0), 1)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1)) + 
  guides(fill = FALSE)

save(p1, file = '~/Downloads/CLL_project/code/figures/CLL_acetylated_pvalue.Rdata')

pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_acetylated_pvalue.pdf'), width = 4, height = 4)
plot(p1)
dev.off()


## Direction plot
aggregated_plots <- lapply(sig_level, function(i) {
  
  countSig <- as.data.frame(resLimma %>% mutate('Dir' = ifelse(logFC > 0, 'Up', 'Down')) %>% group_by(Cluster, Dir, .drop = FALSE) %>% dplyr::filter(adj.P.Val < i) %>% dplyr::count())
})


names(aggregated_plots) <- paste0('FDR - ', sig_level)
aggregated_plots <- reshape2::melt(aggregated_plots)
aggregated_plots$value <- ifelse(aggregated_plots$Dir == 'Up', aggregated_plots$value,
                                 -aggregated_plots$value)
aggregated_plots$Dir <- factor(aggregated_plots$Dir, levels = c('Up', 'Down'))

p <- ggplot(aggregated_plots, aes(x = Cluster, y = value, fill = Dir)) +
  geom_bar(stat = 'identity') +
  scale_y_continuous(labels = abs_axis) +
  scale_fill_manual(name = 'FC direction', values = c('indianred3', 'steelblue3')) +
  facet_grid(.~L1) +
  labs(x = '', y = 'Number of modified differentially abundant peptides', title = 'Acetylation') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  plot.title = element_text(hjust = 0.5))

pdf('~/Downloads/CLL_project/code/figures/CLL_acetylated_FCDirection.pdf', width = 6, height = 5)
plot(p)
dev.off()


# Volcano plot
# Histone gene symbols (HUGO)

resLimmaSig <- resLimma[resLimma$adj.P.Val < levToChoose, ]
toHighlight <- resLimmaSig[grepl(paste(c('^H[0-9]', 'HIST'), collapse = '|'), resLimmaSig$GS), ]

p2 <- ggplot(subset(resLimma, Cluster == 'PG5'), aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(col = 'grey') +
  geom_point(pch = 21, data = subset(resLimmaSig,Cluster == 'PG5'), fill = 'black', size = 2) +
  geom_point(pch = 21, data = subset(toHighlight,Cluster == 'PG5'), fill = 'red', size = 2) +
  geom_vline(xintercept = 0) +
  # geom_vline(xintercept = c(-thresFC, thresFC), lty = 2) +
  geom_hline(yintercept = -log10(levToChoose), lty = 2) +
  # geom_label_repel(data = resLimmaSig, aes(label = GS), min.segment.length = 0.01, size = 3, force = 10) + 
  facet_grid(.~Cluster) +
  scale_y_continuous(expand = c(0,0,0,0.1)) + 
  labs(x = 'Log2 fold change', y = 'Adjusted p-value (-log10)', title = 'Acetylation') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(face = "bold"))


# Change facet colour
p2 <- ggplot_gtable(ggplot_build(p2))
stripr <- which(grepl('strip-t', p2$layout$name))
fills <- as.character(col_list)[5]
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', p2$grobs[[i]]$grobs[[1]]$childrenOrder))
  p2$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf('~/Downloads/CLL_project/code/figures/CLL_acetylated_Volcano_PG5.pdf', width = 4, height = 5)
plot(p2)
dev.off()

pdf('~/Downloads/CLL_project/code/figures/CLL_acetylated_combined.pdf', width = 5, height = 7)
plot.with.inset <-
  ggdraw() +
  draw_plot(p2) +
  draw_plot(p1, x = 0.1, y = 0.5, width = 0.4, height = 0.4)
plot(plot.with.inset)
dev.off()


save(plot.with.inset, file = '~/Downloads/CLL_project/code/figures/CLL_acetylated_combined.RData')



p <- ggplot(resLimma, aes(x = logFC, y = -log10(adj.P.Val))) +
  geom_point(col = 'grey') +
  geom_point(pch = 21, data = resLimmaSig, fill = 'black', size = 2) +
  geom_point(pch = 21, data = toHighlight, fill = 'red', size = 2) +
  geom_vline(xintercept = 0) +
  # geom_vline(xintercept = c(-thresFC, thresFC), lty = 2) +
  geom_hline(yintercept = -log10(levToChoose), lty = 2) +
  # geom_label_repel(data = resLimmaSig, aes(label = GS), min.segment.length = 0.01, size = 3, force = 10) + 
  facet_grid(.~Cluster) +
  scale_y_continuous(expand = c(0,0,0,0.1)) + 
  labs(x = 'Log2 fold change', y = 'Adjusted p-value (-log10)', title = 'Acetylation') + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"), 
        strip.text = element_text(face = "bold"))


# Change facet colour
p <- ggplot_gtable(ggplot_build(p))
stripr <- which(grepl('strip-t', p$layout$name))
fills <- as.character(col_list)
k <- 1
for (i in stripr) {
  j <- which(grepl('rect', p$grobs[[i]]$grobs[[1]]$childrenOrder))
  p$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fills[k]
  k <- k+1
}

pdf('~/Downloads/CLL_project/code/figures/CLL_acetylated_Volcano.pdf', width = 15, height = 4)
plot(p)
dev.off()



## Boxplot for specific peptides
idx <- which(grepl('PG5', resLimmaSig$Cluster) == TRUE)

p_hist <- lapply(idx, function(i) {
  
  box_dat <- data.frame(values = as.numeric(tofit[rownames(tofit) %in% resLimmaSig$Peptide[i], ]),
                        colnames(tofit),
                        cluster = paste0('PG', samples_annot$Cluster))
  
  
  p <- ggplot(box_dat, aes(x = cluster, y = values, fill = cluster)) +
    # geom_violin() +
    geom_boxplot(outlier.colour = NA) +
    geom_jitter(pch = 21, col = 'black', width = 0.1) +
    scale_fill_manual(values = col_list) +
    labs(x = '', y = paste0('Log2 ratio'),
         title = resLimmaSig$GS[i],
         subtitle = resLimmaSig$Gene[i]
    )  +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5, size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 7),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    guides(fill = FALSE)
  
})


pdf('~/Downloads/CLL_project/code/figures/CLL_acetylated_diff_abund_sig.pdf', width = 20, height = 20)
grid.arrange(grobs = p_hist)
dev.off()


## Gene symbol overlap
resLimma$Peptide <- gsub('\\+229.163',  '',  resLimma$Peptide)
resLimma <- resLimma[, c("Peptide", "Gene", "Protein", "GS", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B", "Cluster")]
write.table(resLimma, '~/Downloads/CLL_project/data/CLL_acetylated_diff_abund.txt', sep = '\t', row.names = FALSE)
# write.table(resLimma[resLimma$adj.P.Val < thres, ], '~/Downloads/CLL_project/data/CLL_acetylated_diff_abund_Sig.txt', sep = '\t', row.names = FALSE)


###