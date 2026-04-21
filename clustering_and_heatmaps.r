require(Biobase)
subtype.significant <- read.delim("data/NKI295.subtype.significant.txt", stringsAsFactors=T, row.names=1)
NKI295_sub <- NKI295[rownames(subtype.significant),]
NKI295_sub
NKI295_sub_distM <- dist(t(exprs(NKI295_sub)), method="euclidean")
NKI295_sub_sampleTree <- hclust(NKI295_sub_distM, method="ward.D")
plot(NKI295_sub_sampleTree, cex=0.2)
abline(h=250)
cutree(NKI295_sub_sampleTree, k=4)
cutree(NKI295_sub_sampleTree, h=250)
table(cutree(NKI295_sub_sampleTree, k=4))
table(cutree(NKI295_sub_sampleTree, h=250))
all(names(cutree(NKI295_sub_sampleTree, k=4)) == sampleNames(NKI295_sub))
NKI295_sub$cluster <- cutree(NKI295_sub_sampleTree, k=4)
table(NKI295_sub$cluster)
View(pData(NKI295_sub))
NKI295_sub_sampleTree$order
NKI295_sub$cluster[NKI295_sub_sampleTree$order]
#   install.packages("dendextend")
library(dendextend)
NKI295_sub_sampleTree_dg <- as.dendrogram(NKI295_sub_sampleTree)
subtype.colours <- c("red", "deeppink","dark blue", "light blue","orange")
subtype <- subtype.colours[NKI295_sub$subtype]
ER <- ifelse(NKI295_sub$ER == "Positive", "dark blue", "cornflowerblue")
phenobar <- cbind(subtype, ER)
NKI295_sub_sampleTree_dg <- set(NKI295_sub_sampleTree_dg, "labels", "")
plot(NKI295_sub_sampleTree_dg, axes=F)
colored_bars(colors = phenobar, dend = NKI295_sub_sampleTree_dg, y_shift=-10)
legend("topright", c(levels(NKI295_sub$subtype), paste("ER", levels(NKI295_sub$ER))), fill=c(subtype.colours,"dark blue", "cornflowerblue"), bty="n", cex=1)
NKI295_sub_sampleTree_dg2 <- click_rotate(NKI295_sub_sampleTree_dg)
par(mfrow=c(1,2))
plot(NKI295_sub_sampleTree_dg, axes=F)
colored_bars(colors = phenobar, dend = NKI295_sub_sampleTree_dg,
y_shift=-10)

plot(NKI295_sub_sampleTree_dg2, axes=F)
colored_bars(colors = phenobar, dend = NKI295_sub_sampleTree_dg2, y_shift=-10)
NKI295_sub_distM_t <- dist(exprs(NKI295_sub), method="euclidean")
NKI295_sub_geneTree <- hclust(NKI295_sub_distM_t, method="ward.D")
palette <- colorRampPalette(c("blue", "white", "red"))(32)
heatmap_breaks <- c(-7,seq(-4,4,length.out=31), 7)
heatmap(exprs(NKI295_sub), breaks=heatmap_breaks,
        Colv=as.dendrogram(NKI295_sub_sampleTree),
        Rowv=as.dendrogram(NKI295_sub_geneTree),
        labRow="", labCol="",
        col= palette, scale="none",
        ColSideColors=subtype.colours[NKI295_sub$subtype])
legend("topright",levels(NKI295_sub$subtype), fill= subtype.colours, cex=0.5, bty="n")