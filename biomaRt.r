require(biomaRt)
library("devtools")
devtools::install_version("dbplyr", version = "2.3.4")
mart <- useMart(host="https://feb2014.archive.ensembl.org", biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
ensemblAttributes <- listAttributes(mart)
my_attributes <- c("hgnc_symbol","description","chromosome_name","start_position","end_position","strand","band","ensembl_gene_id") 
my_filters <- "hgnc_symbol" 
my_values = c("ERBB2", "PIK3CA", "BRCA2", "AR", "PTEN")
geneAnnotation <- getBM(attributes=my_attributes, filters=my_filters, values=my_values, mart=mart)
geneAnnotation <- geneAnnotation[!nchar(geneAnnotation$chromosome_name) > 2,]
geneAnnotation <- geneAnnotation[order(as.numeric(geneAnnotation$chromosome_name),geneAnnotation$start_position),]
geneAnnotation$description <- gsub(" $", "", gsub("\\[Sour\\S+\\s+\\S+", "", geneAnnotation$description, perl=T))