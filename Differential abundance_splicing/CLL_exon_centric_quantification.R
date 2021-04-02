### CLL exon centric quant

## Libraries
library(dplyr)
library(limma)
library(ggplot2)


## Data
CLL_metadata <- read.delim('~/Downloads/CLL_project/code/processing/CLL_code/Data/CLL_metadata.txt', sep = '\t', header = TRUE)

# Splicevista output
splicevistaoutput <- read.delim('~/Downloads/CLL_project/code/processing/CLL_code/Data/CLL_splicevista_output.txt', header = TRUE, sep = '\t',  check.names = FALSE)

# Add exon ID
splicevistaoutput$original_exon <- sapply(1:nrow(splicevistaoutput), function(i) {
  paste(unique(c(as.character(splicevistaoutput$start_exon_id[i]),
                 as.character(splicevistaoutput$end_exon_id[i]))), collapse = '_')
})


## Change names
tochange <- grep('POOL', colnames(splicevistaoutput))
set_names <- paste0('Set', rep(c('A', 'B','C', 'D', 'E', 'F', 'G', 'H'), each = 10))
colnames(splicevistaoutput)[tochange] <- gsub('^.*Set[A-Z]_([^ ])','\\1', colnames(splicevistaoutput)[tochange])
colnames(splicevistaoutput)[tochange]  <- paste0(set_names, '_', colnames(splicevistaoutput)[tochange])

# Quant.data
splicevistaoutput_quant <- as.data.frame(splicevistaoutput[, c('Peptide sequence', colnames(splicevistaoutput)[tochange])])
pep_seq <- splicevistaoutput_quant$`Peptide sequence`

total_sample <- ncol(splicevistaoutput_quant) - 1


# Merge replicates
replicates_id <- which(sapply(strsplit(CLL_metadata$TMT, ','), length) == 3)
non_replicates <- splicevistaoutput_quant[, intersect(colnames(splicevistaoutput_quant), CLL_metadata$TMT[-replicates_id])]


replicates <- do.call(cbind ,lapply(replicates_id, function(i) {
  
  rep_id <- trimws(unlist(strsplit(CLL_metadata$TMT[i], ',')))
  dt <- splicevistaoutput_quant[, rep_id]
  dt <- apply(dt, 1, mean, na.rm= TRUE)
  
}))

colnames(replicates) <- CLL_metadata$patient_ID[replicates_id]
colnames(non_replicates) <- CLL_metadata$patient_ID[match(colnames(non_replicates), CLL_metadata$TMT)]

peptide_dt <- cbind.data.frame(non_replicates, replicates)

# peptide_dt$peptide_sequence <- pep_seq
rownames(peptide_dt) <- pep_seq

CLL_metadata <- CLL_metadata[match(colnames(peptide_dt),CLL_metadata$patient_ID), ]
rm(splicevistaoutput_quant)


## Summarize per exon
peptide_dt_exon <- cbind.data.frame('exon' = splicevistaoutput$original_exon, peptide_dt)
exon_dt <- as.data.frame(peptide_dt_exon %>% group_by(exon) %>% summarise_all(median, na.rm = TRUE))
rownames(exon_dt) <- exon_dt$exon
exon_dt$exon <- NULL


# Filter by detection in more that 2 set
exon_dt <- exon_dt[apply(exon_dt, 1, function(i) sum(!is.na(i))) > 18, ]
peptide_dt <- peptide_dt[apply(peptide_dt, 1, function(i) sum(!is.na(i))) > 18, ]
#######################################################################################################
# Limma exon usage
clusters <- CLL_metadata$PG
total_clusters <- length(unique(clusters))


## Differential enrichment analysis
factors <- as.factor(clusters)

design_mat <- model.matrix(~0 + factors)
colnames(design_mat) <- levels(factors)

# Fit linear model
tofit <- exon_dt

fit <- lmFit(tofit, design_mat)
fit <- eBayes(fit, trend=TRUE)

tocompare <- sapply(1:total_clusters, function(i) {
  
  paste0(paste0('PG', i), '-(', paste0('PG', c(1:total_clusters)[-i], collapse = '+'), ')/',total_clusters-1)
})

contrasts <- makeContrasts(contrasts = tocompare, levels =  levels(factors))
fit2 <- contrasts.fit(fit, contrasts)


# Per gene 
geneid <- splicevistaoutput$`Gene(s)`[match(rownames(tofit), splicevistaoutput$original_exon)]
ex <- diffSplice(fit2, geneid=geneid, exonid = rownames(tofit))


exonleveltab <- lapply(1:6, function(i) {
  res <- topSplice(ex, coef = i, number = Inf, test = "t",  sort.by = "none")
  res <- res[complete.cases(res), ]
  res <- res[order(res$FDR), ]
  res$Cluster = paste0('PG', i)
  res
  
})

exonleveltab <- do.call(rbind, exonleveltab)
rownames(exonleveltab) <- NULL

# Thresholds
idSig <- exonleveltab$FDR < 0.01 & abs(exonleveltab$logFC) > 0.5
exonleveltabSig <- exonleveltab[idSig, ]
exonleveltabNonSig <- exonleveltab[!idSig, ]


exonleveltab$GS <- splicevistaoutput$`Associated gene ID(s)`[match(exonleveltab$ExonID,splicevistaoutput$original_exon)]

write.table(exonleveltab[, c("ExonID", "GeneID", "GS", "logFC", "t", "P.Value", "FDR", "Cluster")], '~/Downloads/CLL_project/data/CLL_limma_diffSplice_fdr_001_logFC_05.txt', sep = '\t', row.names = FALSE)

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

## Rename
diffspliceClusterExon$Var1 <- factor(diffspliceClusterExon$Var1, 
                                     levels = c('PG5', 'PG1', 'PG2', 'PG3', 'PG4', 'PG6'), labels = c('ASB-CLL', 'Tris12M-PG', 'Tris12U-PG', 'M-PG', 'U-PG', 'TP53-PG'))

group_cols <- c('Tris12M-PG' = "#ffbe44", 
                'Tris12U-PG' = "#00b35b", 
                'M-PG' = "#9b2193", 
                'U-PG' = "#ec3e6d", 
                'ASB-CLL' = "#9997ff",
                'TP53-PG' = "#753300")

p <- ggplot(diffspliceClusterExon, aes(x = Var1, y = Freq, fill = Var1)) + 
  geom_bar(stat = 'identity') + 
  geom_text(mapping = aes(label  = fisher), size = 5) + 
  scale_fill_manual(values = group_cols) +
  labs(x = '', y = 'Nr. exons with differential exon usage \n(Adjusted p < 0.01 and |log2FC| > 0.5)', title = 'Peptides - exon usage') +  
  # subtitle = 'Exon level') + 
  scale_y_continuous(expand = c(0,0, 0.1, 0)) + 
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5), 
        plot.subtitle = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45, hjust = 1), 
        panel.grid = element_blank()) + 
  guides(fill = FALSE)

# save(p, file = '~/Downloads/CLL_project/code/figures/CLL_diffspl_Exon_usage_peptides.Rdata')
# pdf(paste0('~/Downloads/CLL_project/code/figures/CLL_diffSplice_barplot_exon.pdf'), width = 4, height = 6)
plot(p)
# dev.off()

###