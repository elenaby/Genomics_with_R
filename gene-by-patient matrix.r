TCGA_id <- "GBM"
mutations <- read.delim(paste0("data/TCGA/",TCGA_id , ".mutations.txt"), stringsAsFactors=F)
mutations <- mutations[mutations$Variant_Classification != "Silent",]
mutations <- mutations[mutations$Hugo_Symbol != "Unknown",]
mutation_matrix <- matrix(NA, nrow=length(unique(mutations$Hugo_Symbol)), ncol=length(unique(mutations$Tumor_Sample_Barcode)))
colnames(mutation_matrix) <- unique(mutations$Tumor_Sample_Barcode)
row.names(mutation_matrix) <- unique(mutations$Hugo_Symbol)
head(mutation_matrix[,1:3])
for(tumourID in colnames(mutation_matrix)){
    mutation_matrix[,tumourID] <- row.names(mutation_matrix) %in% mutations$Hugo_Symbol[mutations$Tumor_Sample_Barcode == tumourID]
}
mutation_matrix <- mutation_matrix[order(rowSums(mutation_matrix), decreasing=T),]
mutation_matrix <- mutation_matrix[row.names(mutation_matrix) != "Unknown",]
barplot(100*apply(mutation_matrix[1:30,], 1, sum)/ncol(mutation_matrix), las=2, cex.names=0.75, main=paste(TCGA_id, "mutations"), ylab="frequency", col="dark green")
write.table(mutation_matrix, paste0("results/",TCGA_id , ".mutation_matrix.matrix.xls"), row.names=F, sep="\t", na="")
