BiocManager::install("fgsea")
library(fgsea)
HER2_ranks <- brca_HER2_res$log2FoldChange
names(HER2_ranks) <- brca_HER2_res$symbol
HER2_ranks <- HER2_ranks[!is.na(HER2_ranks)]
pathways <- gmtPathways("data/gene_sets/c2.all.v7.0.symbols.gmt")
fgseabrca_HER2_res <- fgsea(pathways, HER2_ranks, nperm=10000, maxSize=500, minSize=15)
fgseabrca_HER2_res <- fgseabrca_HER2_res[fgseabrca_HER2_res$pval < 0.05,]
fgseabrca_HER2_res <- fgseabrca_HER2_res[order(fgseabrca_HER2_res$NES, decreasing=T),]
head(fgseabrca_HER2_res[,1:7], 10)
library(ggplot2)
plotEnrichment(pathways[[fgseabrca_HER2_res$pathway[1]]],HER2_ranks) + ggtitle(fgseabrca_HER2_res$pathway[1])
par(mar=c(5,20,2,2))
barplot(fgseabrca_HER2_res$NES[20:1], horiz=T, names.arg=fgseabrca_HER2_res$pathway[20:1], las=1, cex.names=0.6, xlim=c(2,3.5), xpd=F, xlab="NES", col="dark blue")
title("fgsea \nBRCA HER2 enrichments")
fgseabrca_HER2_resdf <- as.data.frame(fgseabrca_HER2_res[,1:7])
fgseabrca_HER2_resdf $leadingEdge <- unlist(lapply(fgseabrca_HER2_res $leadingEdge, function(x) paste(x, collapse=",")))

write.table(fgseabrca_HER2_resdf, "results/fgsea_brca_HER2_results.xls", sep="\t", na="", row.names=F, quote=F)