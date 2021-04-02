### CLL trisomy12/IGHV status interaction

## Libraries
library(biomaRt)
library(DEqMS)
library(ggplot2)
library(ggrepel)
library(GSA)
library(fgsea)

## Data
msig <- GSA.read.gmt('Data/c2.cp.kegg.v6.2.symbols.gmt.txt')
enrich_list <- msig$genesets
names(enrich_list) <- msig$geneset.names


# Colors
col_list <- c('PG1' = "#ffbe44",'PG2'="#00b35b",'PG3'="#9b2193",'PG4' = "#ec3e6d",
              'PG5' = "#9997ff",'PG6' = "#753300")

CLL_metadata <- read.delim('Data/CLL_metadata.txt', sep = '\t', row.names = 1)

psm_per_set <- read.delim('Data/symbols_table.txt', sep = '\t')
psm_per_set_quant <- psm_per_set[, grep('quanted_psm_count', colnames(psm_per_set))]
rownames(psm_per_set_quant) <- psm_per_set$Gene.Name

# Minimun across sets
psm_per_set_quant[is.na(psm_per_set_quant)] <- 0
psm_per_set_quant$count <- apply(psm_per_set_quant, 1, min, na.rm = TRUE)
psm_per_set_quant <- psm_per_set_quant[psm_per_set_quant$count != 0, ]


load('Data/CLL_multiomics_MAE_19022021.RData')
proteomics <- multiomics_MAE@ExperimentList$proteomics
transcriptomics <- multiomics_MAE@ExperimentList$RNAseq_norm


# Convert ensemb to gene symbols
ensembl_92 <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                      host="http://apr2018.archive.ensembl.org", 
                      path="/biomart/martservice",
                      dataset="hsapiens_gene_ensembl")

annot <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'chromosome_name'),
               mart = ensembl_92, 
               filters = 'hgnc_symbol', 
               values = rownames(proteomics), 
               useCache = FALSE)

# Remove NA
annot <- annot[annot$hgnc_symbol != '', ]
annot <- annot[!duplicated(annot$ensembl_gene_id), ]

proteomics <- proteomics[rownames(proteomics) %in% annot$hgnc_symbol, ]
rownames(proteomics) <- annot$ensembl_gene_id[match(rownames(proteomics), annot$hgnc_symbol)]

psm_per_set_quant <- psm_per_set_quant[rownames(psm_per_set_quant) %in% annot$hgnc_symbol, ]
rownames(psm_per_set_quant) <- annot$ensembl_gene_id[match(rownames(psm_per_set_quant), annot$hgnc_symbol)]

common_genes <- Reduce(intersect, list(rownames(transcriptomics), 
                                       rownames(proteomics), 
                                       rownames(psm_per_set_quant)))

common_samples <- Reduce(intersect, list(
  colnames(transcriptomics),
  colnames(proteomics), 
  rownames(CLL_metadata)))


## Spearman correlations per Trisomy group
trisomy_12_pos <- common_samples[which(CLL_metadata[common_samples, 'trisomy12'] == 1)]
trisomy_12_neg <- common_samples[which(CLL_metadata[common_samples, 'trisomy12'] == 0)]



## Subset samples
tokeep <- !is.na(CLL_metadata$IGHV) & !is.na(CLL_metadata$trisomy12)
CLL_metadata_interaction <- CLL_metadata[tokeep, ]


proteomics_samples <- intersect(rownames(CLL_metadata_interaction), colnames(proteomics))
CLL_metadata_interaction <- CLL_metadata_interaction[proteomics_samples, ]
proteomics_interaction <- proteomics[,proteomics_samples]
#####################################################################################
# Linear regression
regression_samples <- intersect(common_samples, proteomics_samples)

CLL_metadata_regress  <- CLL_metadata[regression_samples, ]

# Scale
scaled_proteomics <- as.data.frame(t(scale(t(proteomics[, regression_samples]))))
scaled_transcriptomics <- as.data.frame(t(scale(t(transcriptomics[, regression_samples]))))

rownames(scaled_proteomics) <- rownames(proteomics)
rownames(scaled_transcriptomics) <- rownames(transcriptomics)

coef_lm <- sapply(common_genes, function(j) {
  
  res_lm <- lm(as.numeric(scaled_proteomics[j,]) ~  as.numeric(scaled_transcriptomics[j,]) + CLL_metadata_regress$IGHV*CLL_metadata_regress$trisomy12)
  res_lm <- summary(res_lm)
  res_lm$coefficients[2,1]
  })

coef_dt <- data.frame('ENS_ID' = common_genes, 'Symbol' = annot$hgnc_symbol[match(common_genes, annot$ensembl_gene_id)], 'coef' = coef_lm)


coef_dt <- coef_dt[order(-coef_dt$coef), ]
# write.csv(coef_dt, 'Data/CLL_IGHV_Tri12_beta.csv',row.names = FALSE)
#####################################################################################
## Differential abundance 
CLL_metadata_interaction$IGHV <- factor( CLL_metadata_interaction$IGHV , levels = c('U', 'M'), labels = c('IGHV-U', 'IGHV-M'))
# CLL_meta_sub$Group5 <- ifelse(CLL_meta_sub$CCP2 == 5, 'Group5', 'Rest')
CLL_metadata_interaction$trisomy12 <- factor( CLL_metadata_interaction$trisomy12, levels = c(0,1), labels = c('Trisomy 12-', 'Trisomy 12+'))

## DEqMS / Limma test
design_mat <- model.matrix(~CLL_metadata_interaction$IGHV + CLL_metadata_interaction$trisomy12 + CLL_metadata_interaction$IGHV:CLL_metadata_interaction$trisomy12)

# DEqMS - protein
tofit <- proteomics_interaction
fit <- lmFit(tofit, design_mat)
fit2 <- eBayes(fit)

fit2$count = psm_per_set_quant[rownames(proteomics_interaction), 'count']
fit3 <- spectraCounteBayes(fit2)

# Extact adjusted pvalues
resDeqMS = outputResult(fit = fit3, coef_col = 4)
resDeqMS$eff_size <- fit3$coefficients[rownames(resDeqMS), 4]/sqrt(fit3$sca.postvar)[rownames(resDeqMS)]
# resDeqMS$eff_size <- scale(resDeqMS$logFC)
resDeqMS$Gene <- rownames(resDeqMS)
resDeqMS$Symbol <- annot$hgnc_symbol[match(rownames(resDeqMS), annot$ensembl_gene_id)]
resDeqMS$chromosome <- annot$chromosome_name[match(rownames(resDeqMS), annot$ensembl_gene_id)]
################################################################################################################
thresP <- 0.05
thresFC <- 1

## Volcano plot
tohighlight <- resDeqMS[abs(resDeqMS$logFC) > thresFC & resDeqMS$adj.P.Val < thresP, ]
tohighlight_MAPK <- resDeqMS[resDeqMS$Symbol %in% enrich_list[['KEGG_B_CELL_RECEPTOR_SIGNALING_PATHWAY']], ]

p <- ggplot(resDeqMS, aes(x = logFC, y = -log10(adj.P.Val), size = abs(logFC) *  -log10(adj.P.Val))) +
  geom_vline(xintercept = c(0, -thresFC, thresFC), lty = c(1,2,2)) + 
  geom_hline(yintercept = -log10(thresP), lty = 2) +
  annotate(geom = 'label', fill = 'red', col ='white', x = c(-thresFC, thresFC), y = 3.5, 
           label = c('Higher in IGHV-U', 'Higher in IGHV-M')) +
  geom_point() + 
  geom_point(pch  = 21, data = tohighlight, fill = 'red', col = 'black')  + 
  geom_point(pch  = 21, data = tohighlight_MAPK, fill = 'yellow', col = 'black')  + 
  geom_text_repel(tohighlight, mapping = aes(label = Symbol, size = abs(logFC) *  -log10(adj.P.Val)/1.5), min.segment.length = 0.01,col = 'red') + 
 
  scale_y_continuous(expand =  c(0,0, 0.1,0)) + 
  theme_bw() + 
  labs(x = 'log2 fold change', y = '-log10 adjusted p-value', title = 'IGHV:Trisomy12 interaction - Proteomics') + 
  theme(plot.title = element_text(hjust = 0.5), 
        panel.grid = element_blank()) + 
  guides(size = FALSE)

pdf(file = paste0("Figures/CLL_IGHV_Tri12_volcano.pdf"), width = 4, height = 6)
plot(p)
dev.off()

deqms_tb <- resDeqMS[, c("Gene", "Symbol", "chromosome" , "logFC", "AveExpr", "t", "P.Value"    ,"adj.P.Val",   "B",  "count",  "sca.t", "sca.P.Value", "sca.adj.pval", "eff_size")]
deqms_tb <- deqms_tb[order(-deqms_tb$eff_size), ]

# write.csv(deqms_tb, file = 'Data/CLL_deqMS_IGHV_tri12.csv', row.names = FALSE)
#################################################################################
## Enrichment analysis
genesToCheck <- resDeqMS$logFC
names(genesToCheck) <- resDeqMS$Symbol
genesToCheck <- sort(genesToCheck,  decreasing = TRUE)

# Over-representation analysis
set.seed(1234)
enrich_dt <- fgseaMultilevel(pathways = enrich_list,
                    stats = genesToCheck,
                    minSize=15,
                    maxSize=200)

enrich_dt <- as.data.frame(enrich_dt)


# ## T-test of leading edges
# tokeep <- sapply(1:nrow(enrich_dt), function(i) {
#   leading_edge <- enrich_dt$leadingEdge[[i]]
#   leading_edge <- intersect(leading_edge, coef_dt$Symbol)
#   length(enrich_dt$leadingEdge[[i]]) > 4
# })
# enrich_dt <- enrich_dt[tokeep, ]

enrich_dt[, c(
  # 'cohens_d', 
  't_test_pval', 'median_coef', 'n')] <- t(sapply(1:nrow(enrich_dt), function(i) {
  # print(i)
  leading_edge <- enrich_dt$leadingEdge[[i]]
  leading_edge <- intersect(leading_edge, coef_dt$Symbol)
  if(length(leading_edge) < 5) {
    return(rep(NA, 3))
  }
  

  t_res <- t.test(coef_dt$coef[coef_dt$Symbol %in% leading_edge],
                  coef_dt$coef[!coef_dt$Symbol %in% leading_edge])
  # cohens_d <- cohen.d(coef_dt$coef[coef_dt$Symbol %in% leading_edge], coef_dt$coef[!coef_dt$Symbol %in% leading_edge])
  c(
    # cohens_d$estimate, 
    t_res$p.value, 
    median(coef_dt$coef[coef_dt$Symbol %in% leading_edge]),
    length(leading_edge))
  
}))

enrich_dt$t_test_adjpval <- p.adjust(enrich_dt$t_test_pval, method = 'BH')
enrich_dt$cor_dir <- ifelse(enrich_dt$median_coef < 0.23 & enrich_dt$t_test_adjpval < 0.01, 'Down', ifelse(enrich_dt$median_coef > 0.23 & enrich_dt$t_test_adjpval < 0.01, 'Up', ifelse(enrich_dt$t_test_adjpval > 0.01, 'Not sig.', NA)))

enrich_dt$cor_dir <- factor(enrich_dt$cor_dir, levels = c('Up', 'Down', 'Not sig.'))


p <- ggplot(enrich_dt, aes(x = NES, y = -log10(padj), col = cor_dir, size = -log10(padj) * abs(NES),  fill =  -log10(padj))) + 
  geom_point(pch = 21, stroke = 2) + 
  geom_text_repel(data = subset(enrich_dt, padj < 0.01 & cor_dir %in% c('Up' ,'Down')), mapping = aes(label = pathway), size =2, show.legend = FALSE) + 
  scale_fill_distiller(palette = "YlOrRd") + 
  scale_color_manual(name = 'Protein-RNA corr.', 
                     values = c('Up' = 'indianred3', 'Down' = 'steelblue3', 'Not sig.' = 'grey'), na.value= 'black') + 
  theme_bw() + 
  theme(plot.background = element_blank()) +
  guides(size = FALSE, fill = FALSE) 



pdf(file = paste0("Figures/CLL_IGHV_Tri12_volcano_enrich.pdf"), width = 8, height = 9)
plot(p)
dev.off()

## Save data
# enrich_dt <- enrich_dt[, c("pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge", "n", 'cor_dir', "t_test_pval", "t_test_adjpval")]
# 
# 
# enrich_dt$leadingEdge <- sapply(enrich_dt$leadingEdge, paste, collapse = ',')
# enrich_dt <- enrich_dt[order(enrich_dt$padj), ]
# 
# 
# write.csv(enrich_dt, file = 'Data/CLL_IGHV_Tri12_enrichment_results.csv', row.names = FALSE)

###