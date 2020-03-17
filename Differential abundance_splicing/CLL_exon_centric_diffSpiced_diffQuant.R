### CLL exon centric differential splicing

## Libraries
# install.packages("easypackages")
library(easypackages)

toInstall <- c('limma', 'ggplot2', 'reshape2', 'dplyr')

# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install(toInstall)

libraries(toInstall)


# Limma splice variants
# SetIDs
samples_annot <- read.table('~/Downloads/CLL_project/data/patIDmeta_translation_table.txt', header = TRUE)

# Correct group
meta <- read.table('~/Downloads/CLL_project/data/CLL_Meta_Proteomics_2019-12-09.txt', sep = '\t', header = TRUE)

# Change according to the new cluster
samples_annot$Cluster[match(meta$Proteomics_Pat_ID, samples_annot$Labels_proteomics)] = meta$CCP2


# Data
# Splicevista output
CLL_exon <- read.table('~/Downloads/CLL_project/data/CLL_exon_added_transcripts.txt', header = TRUE, sep = '\t',  check.names = FALSE)

# Merged replicates by average,  grouped exons by peptide median (exonic level)
non_replicates_samples <- read.table('~/Downloads/CLL_project/data/CLL_exon_samples.txt', sep = '\t', header = TRUE, check.names = FALSE)

# Peptide sequence level
non_replicates <- read.table('~/Downloads/CLL_project/data/CLL_sequence_samples.txt', sep = '\t',  header = TRUE, check.names = FALSE)


## Filter by detection in more that 2 set
non_replicates_samples <- non_replicates_samples[apply(non_replicates_samples, 1, function(i) sum(!is.na(i))) > 18, ]
non_replicates <- non_replicates[apply(non_replicates[,1:68], 1, function(i) sum(!is.na(i))) > 18, ]
non_replicates$Transcript <- as.character(non_replicates$Transcript)


# Pairwise differential abundance
clusters <- paste0('PG', meta$CCP2)


## Differential enrichment analysis
factors <- as.factor(clusters)

design_mat <- model.matrix(~0 + factors)
colnames(design_mat) <- levels(factors)

# Fit linear model
tofit <- non_replicates_samples
fit <- lmFit(tofit, design_mat)

tocompare <- c("PG1-(PG2+PG3+PG4+PG5+PG6)/5",
               "PG2-(PG1+PG3+PG4+PG5+PG6)/5", 
               "PG3-(PG1+PG2+PG4+PG5+PG6)/5", 
               "PG4-(PG1+PG2+PG3+PG5+PG6)/5",
               "PG5-(PG1+PG2+PG3+PG4+PG6)/5", 
               "PG6-(PG1+PG2+PG3+PG4+PG5)/5")


clusters <-  levels(factors)
contrasts <- makeContrasts(contrasts = tocompare, levels =  clusters)

fit2 <- contrasts.fit(fit, contrasts)


geneid <- CLL_exon$`Gene(s)`[match(rownames(tofit), paste(CLL_exon$start_exon_id, CLL_exon$end_exon_id, sep = '_'))]
ex <- diffSplice(fit2, geneid=geneid, exonid = rownames(non_replicates_samples))

## Per gene 

# # Barplot of differential exon usage
# res_gene <- lapply(1:6, function(i) {
#   resTop <- topSplice(fit = ex, coef = i, number=nrow(tofit), sort.by = "p")
#   resTop <- resTop[complete.cases(resTop), ]
#   
#   resTop$Cluster <- paste0('PG', i)
#   
#   resTop
# })
# 
# res_gene <- do.call(rbind, res_gene)
# 
# 
# # Fisher's exact test
# countSig <- as.data.frame(res_gene %>% group_by(Cluster, .drop=FALSE) %>% dplyr::filter(FDR < 0.01 ) %>% dplyr::count())
# countNonSig <- as.data.frame(res_gene %>% group_by(Cluster, .drop=FALSE) %>% dplyr::filter(FDR > 0.01) %>% dplyr::count())
#   
# resp <- sapply(1:6, function(j) {
#   tab <- matrix(c(countSig$n[j],sum(countSig$n) - countSig$n[j], 
#                     countNonSig$n[j],sum(countNonSig$n) - countNonSig$n[j]),
#                   byrow = TRUE,
#                   nrow = 2,
#                   dimnames = list(FC = c("Mod", "NonMod"),
#                                   Group = c('This', 'Rest')))
#     res <- fisher.test(tab, alternative = 'greater')
#     res$p.value
#     
#   })
# resp <- p.adjust(resp, method = 'BH')
# fisherTab <- data.frame('Sig' = ifelse(resp< 0.001, '***', ifelse(resp< 0.01, '**', ifelse(resp< 0.1, '*', ''))), 
#                         'n' = countSig$n,
#                         'Cluster' = paste0('PG',1:6))
#   
# 
# p1 <- ggplot(countSig, aes(x = Cluster, y = n, fill = Cluster)) + 
#   geom_bar(stat = 'identity') + 
#   geom_text(fisherTab, mapping = aes(label = Sig), size = 5) + 
#   scale_fill_manual(values  = c("#ffbe44", "#00b35b", "#9b2193", "#ec3e6d", "#9997ff", "#753300")) +
#   labs(x = '', y = 'Number of genes with differential exon usage (Adjusted p-value < 0.01)', title = 'DiffSplice', subtitle = 'Gene level') +
#   scale_y_continuous(expand = c(0, 0, 0.1, 0)) + 
#   theme_bw() + 
#   theme(panel.grid = element_blank(),
#         strip.background = element_rect(fill = "white"), 
#         strip.text = element_text(face = "bold"),
#         plot.title = element_text(hjust = 0.5), 
#         plot.subtitle = element_text(hjust = 0.5),
#         axis.text.x = element_text(angle = 45, hjust = 1)) + 
#   guides(fill = FALSE)
# 
# save(p1, file = '~/Downloads/CLL_project/code/figures/CLL_diffspl_Gene.Rdata')


# pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_diffSplice_barplot_gene.pdf'), width = 4, height = 6)
# plot(p)
# dev.off()

# Plot Exon level 
exonleveltab <- lapply(1:6, function(i) {
  res <- topSplice(ex, coef = i, number = Inf, test = "t",  sort.by = "none")
  res <- res[complete.cases(res), ]
  res <- res[order(res$FDR), ]
  res$Cluster = paste0('PG', i)
  res
  
  })

exonleveltab <- do.call(rbind, exonleveltab)
rownames(exonleveltab) <- NULL

idSig <- exonleveltab$FDR < 0.01 & abs(exonleveltab$logFC) > 0.5
exonleveltabSig <- exonleveltab[idSig, ]
exonleveltabNonSig <- exonleveltab[!idSig, ]


exonleveltab$GS <- CLL_exon$`Associated gene ID(s)`[match(exonleveltab$ExonID, paste(CLL_exon$start_exon_id, CLL_exon$end_exon_id, sep = '_'))]

write.table(exonleveltab[, c("ExonID", "GeneID", "GS", "logFC", "t", "P.Value", "FDR", "Cluster")], '~/Downloads/CLL_project/data/CLL_limma_diffSplice_fdr_001.txt', sep = '\t', row.names = FALSE)

diffspliceClusterExonNon <-  as.data.frame(table(exonleveltabNonSig$Cluster))
diffspliceClusterExon <- as.data.frame(table(exonleveltabSig$Cluster))
diffspliceClusterExon$Var1 <- gsub('Cluster', '', diffspliceClusterExon$Var1)

diffspliceClusterExon$fisher <- sapply(1:6, function(j) {
  
  tab <- matrix(c(diffspliceClusterExon$Freq[j],sum(diffspliceClusterExon$Freq) - diffspliceClusterExon$Freq[j], 
                  diffspliceClusterExonNon$Freq[j],sum(diffspliceClusterExonNon$Freq) - diffspliceClusterExonNon$Freq[j]),
                byrow = TRUE,
                nrow = 2,
                dimnames = list(FC = c("Sig", "NonSig"),
                                Group = c('This', 'Rest')))
  res <- fisher.test(tab, alternative = 'greater')
  res$p.value
  
})
diffspliceClusterExon$fisher <- p.adjust(diffspliceClusterExon$fisher, method = 'BH')
diffspliceClusterExon$fisher <- ifelse(diffspliceClusterExon$fisher < 0.001, '***', ifelse(diffspliceClusterExon$fisher < 0.01, '**', ifelse(diffspliceClusterExon$fisher < 0.1, '*', '')))


p2 <- ggplot(diffspliceClusterExon, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = 'identity') + 
  geom_text(mapping = aes(label  = fisher), size = 5) + 
  scale_fill_manual(values  = c("#ffbe44", "#00b35b", "#9b2193", "#ec3e6d", "#9997ff", "#753300")) +
  labs(x = 'PG', y = 'Number of exons with differential exon usage \n(Adjusted p-value < 0.01 and |log2FC| > 0.5)', title = 'DiffSplice', subtitle = 'Exon level') + scale_y_continuous(expand = c(0,0, 0.1, 0)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0), 
        panel.grid = element_blank()) + 
  guides(fill = FALSE)

save(p2, file = '~/Downloads/CLL_project/code/figures/CLL_diffspl_Exon.Rdata')


pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_diffSplice_barplot_exon.pdf'), width = 4, height = 6)
plot(p2)
dev.off()


## Normalize by differential abundance
# install.packages("easypackages")
fit3 <- eBayes(fit2)

resLimma = lapply(1:length(clusters), function(i) {
  res = topTable(fit3, coef = i, number = Inf)
  res <- data.frame('Mapped_transcripts' = rownames(res), 'FC' = res$logFC, 'AdjPval' = res$adj.P.Val, 'Cluster' = paste0('PG', i))
  
  # Add meta
  res <- cbind.data.frame(res, CLL_exon[match(res$Mapped_transcripts, paste(CLL_exon$start_exon_id, CLL_exon$end_exon_id, sep = '_')), c("Protein(s)", "Gene(s)", "Associated gene ID(s)", "chr", "start", "end", "strand", "pep_start_exon", "pep_end_exon", "start_exon_id", "end_exon_id", "NumOf.coding_transcript.From.gene", "NumOf.trans.mapped.to", "Peptide sequence", 'Transcript')])
  
  # Sort by pvalue
  res <- res[order(res$AdjPval), ]
  res <- res[res$AdjPval < 0.01 & abs(res$FC) > 0.5, ]
  res <- res[complete.cases(res), ]
  
  res
})

resLimma <- do.call(rbind, resLimma)
rownames(resLimma) <- NULL


## Per cluster siginificantly differentiated exons 
diffabundfreq <- as.data.frame(table(resLimma$Cluster))

p3 <- ggplot(diffabundfreq, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = 'identity') + 
  scale_fill_manual(values  = c("#ffbe44", "#00b35b", "#9b2193", "#ec3e6d", "#9997ff", "#753300")) +
  labs(x = '', y = 'Number of exons with differential abundance \n(Adjusted p-value < 0.01 and |log2FC| > 0.5)', title = 'Limma',
       subtitle = 'Exon level') +
  scale_y_continuous(expand = c(0,0, 0.1, 0)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank()) + 
  guides(fill = FALSE)

save(p3, file = '~/Downloads/CLL_project/code/figures/CLL_diffabund_diffspl.Rdata')


# pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_diffLimma_barplot_exon.pdf'), width = 4, height = 6)
# plot(p3)
# dev.off()

perClusterOverlap <- sapply(1:6, function(i) {
  spl <- unique(exonleveltabSig$GeneID[exonleveltabSig$Cluster == paste0('PG', i)])
  abun <- unique(resLimma$`Gene(s)`[resLimma$Cluster == paste0('PG', i)])
  dat = c('abundant' = length(abun), 'abundant_spliced' = length(intersect(spl, abun)))
  # list(id = intersect(spl, abun), len = length(intersect(spl, abun)), per = length(intersect(spl, abun))/length(abun))
  
})


# Normalized by differentially abundant
tab <- unlist(perClusterOverlap['abundant_spliced', ]/perClusterOverlap['abundant', ])

resp <- sapply(1:6, function(j) {
  tab <- matrix(c(perClusterOverlap['abundant_spliced',j],sum(perClusterOverlap['abundant_spliced',]) - perClusterOverlap['abundant_spliced',j], 
                  perClusterOverlap['abundant',j], sum(perClusterOverlap['abundant', ]) - perClusterOverlap['abundant',j]),
                byrow = TRUE,
                nrow = 2,
                dimnames = list(FC = c("Spliced", "Abundant"),
                                Group = c('This', 'Rest')))
  res <- fisher.test(tab, alternative = 'greater')
  res$p.value
  
})
resp <- p.adjust(resp, method = 'BH')


tab <- melt(tab)
tab$Cluster <-  paste0('PG', 1:6)

fisherTab <- data.frame('Sig' = ifelse(resp< 0.001, '***', ifelse(resp< 0.01, '**', ifelse(resp< 0.1, '*', ''))), 
                        'value' = tab$value,
                        'Cluster' = paste0('PG',1:6))

p4 <- ggplot(tab, aes(x = Cluster, y = value, fill = Cluster)) + 
  geom_bar(stat = 'identity') + 
  geom_text(aes(label = round(value, 2)), nudge_y = 0.05) + 
  geom_text(fisherTab, mapping = aes(label = Sig), size = 5, hjust = -0.1) + 
  scale_fill_manual(values  = c("#ffbe44", "#00b35b", "#9b2193", "#ec3e6d", "#9997ff", "#753300")) +
  labs(title = 'Differential exon abundance-usage normalized by differential exon abundance',
       subtitle = 'Gene level') + 
  scale_y_continuous(expand = c(0,0, 0.1, 0)) + 
  coord_polar() + 
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5), 
        panel.grid = element_blank(),
        strip.background = element_rect(fill = 'white')) + 
  guides(fill = FALSE)


save(p4, file = '~/Downloads/CLL_project/code/figures/CLL_diffabund_diffspl_Polar.Rdata')


# pdf('~/Downloads/CLL_project/code/figures/Exon_quantification_circle_plot.pdf', width = 6, height = 6)
# plot(p)
# dev.off()
####################################################################################################