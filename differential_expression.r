library(DESeq2)
brca_counts <- read.delim("data/tcga_brca_counts_table.txt", stringsAsFactors=F )
tcga_brca_clin <- read.delim("data/tcga_brca_clinical_RNASeq.txt",
stringsAsFactor=T, na.strings=c("", "NA", " "))
brca_HER2_dds <- DESeqDataSetFromMatrix(
    countData = brca_counts,
 colData = tcga_brca_clin,
 design= ~ HER2_STATUS)
brca_HER2_dds <- DESeq(brca_HER2_dds)
brca_HER2_res <- as.data.frame(results(brca_HER2_dds))
brca_HER2_res <- brca_HER2_res[order(brca_HER2_res$padj),]
library(biomaRt)
mart <- useMart(host="feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
annotation <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol", "chromosome_name", "start_position", "end_position"), filters="ensembl_gene_id", values = row.names(brca_HER2_res), mart=mart)
names(annotation) <- c("ENSG", "symbol", "chrom", "start", "end")
brca_HER2_res <- merge(annotation, brca_HER2_res, by="ENSG",all.x=F, all.y=T)
brca_HER2_res <- brca_HER2_res[!is.na(brca_HER2_res$symbol),]
brca_HER2_res <- brca_HER2_res[brca_HER2_res$symbol != "",]
brca_HER2_res <- brca_HER2_res[order(brca_HER2_res$padj),]
brca_HER2_res <- brca_HER2_res[!duplicated(brca_HER2_res$symbol),
View(brca_HER2_res)]