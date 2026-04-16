GBM.vcf <- read.delim("data/GBM.vcf", header=T, stringsAsFactors=F)
head(GBM.vcf, 2)
ensembl71 <- read.delim("data/ensembl71.txt", sep="\t", stringsAsFactors=F, header=T)
GBM.vcf <- merge(ensembl71, GBM.vcf, by.x="ensg", by.y="Gene", all.x=T, all.y=F)