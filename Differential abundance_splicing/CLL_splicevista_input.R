### CLL exon level identification - create input files

## Libraries
library(data.table)
library(biomaRt)
library(seqinr)


## BiomaRt (version 92)
ensembl <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                  host="http://apr2018.archive.ensembl.org", 
                  path="/biomart/martservice",
                  dataset="hsapiens_gene_ensembl")

annot <- getBM(attributes=c('ensembl_gene_id', 'ensembl_transcript_id', 'ensembl_peptide_id_version'),
mart = ensembl)


# Remove unmapped peptides
annot <- annot[annot$ensembl_peptide_id_version != '', ]


write.table(x = annot, file = 'Data/martR_GRCh38_version.txt', sep = '\t', row.names = FALSE, quote = FALSE)


## Peptide table
peptide_table <- fread('Data/Proteomics - search output/HiRIEF/peptides_table.txt', sep = '\t')
toindex <- colnames(peptide_table)[grep('POOL', colnames(peptide_table))]
peptide_table <- peptide_table[, c("Peptide sequence", "Protein(s)", "Gene(s)", "Associated gene ID(s)", toindex), with=FALSE]

peptide_table$`Peptide sequence` <- gsub("[^A-Z]", "", peptide_table$`Peptide sequence`)

# Overlap with biomart 
peptide_table <- peptide_table[peptide_table$`Protein(s)` %in% annot$ensembl_peptide_id_version & peptide_table$`Gene(s)` %in% annot$ensembl_gene_id, ]


# Merge peptides by median
peptide_table <- peptide_table[,lapply(.SD, median, na.rm = TRUE), by = c('Peptide sequence', 'Protein(s)', "Gene(s)", "Associated gene ID(s)")]


peptide_table <- as.data.frame(peptide_table)


# Remove complete NAs
na_num <- apply(peptide_table[,-c(1:4)], 1, function(i) sum(is.na(i)))
peptide_table <- peptide_table[-which(na_num == (ncol(peptide_table) - 4)), ]


write.table(x = peptide_table, file = 'Data/CLL_pep_collapsed.txt', sep = '\t', quote = FALSE, row.names = FALSE)


# Splicevista command - Generate SpliceVista output file
# python SpliceVista.py --input CLL_pep_collapsed.txt  --gtf Homo_sapiens.GRCh38.92.gtf --fasta Homo_sapiens.GRCh38.pep.all.fa  --IDmap martR_GRCh38_version.txt --output CLL_output.txt
###