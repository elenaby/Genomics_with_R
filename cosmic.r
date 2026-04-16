cosmic <- read.delim("data/CosmicMutantExportCensus.tsv", stringsAsFactors=F)
cosmic.PIK3CA <- cosmic[cosmic$Gene.name == "PIK3CA",]
write.table(cosmic.PIK3CA, "results/COSMIC.PIK3CA.xls", sep="\t", na="", row.names=F)
PIK3CA.mutations <- data.frame(table(cosmic.PIK3CA$Mutation.AA))
names(PIK3CA.mutations) <- c("Mutation.AA", "Count")
PIK3CA.mutations$Frequency <- 100*PIK3CA.mutations$Count/sum(PIK3CA.mutations$Count)
PIK3CA.mutations  <- PIK3CA.mutations[order(PIK3CA.mutations$Count, decreasing=T),]
write.table(PIK3CA.mutations, "results/PIK3CA.COSMIC.mutation.counts.xls", sep="\t", row.names=F, na="", quote=F)
PIK3CA.mutations$protein.position <- as.numeric(
    gsub("p\\..([0-9]+).+", "\\1", as.character(PIK3CA.mutations$Mutation.AA))
)
PIK3CA.mutations <- PIK3CA.mutations[order(PIK3CA.mutations$Frequency, decreasing=T),]
plot(PIK3CA.mutations$protein.position, PIK3CA.mutations$Frequency, typ="h")
plot(PIK3CA.mutations$protein.position, PIK3CA.mutations$Frequency, typ="h", lwd=2, main="PIK3CA", xlab="Protein Position", ylab="Frequency", ylim=c(0,max(PIK3CA.mutations$Frequency)*1.5))

text(PIK3CA.mutations$protein.position[1:3], PIK3CA.mutations$Frequency[1:3], labels=PIK3CA.mutations$Mutation.AA[1:3], srt=90, adj=c(-0.1,0.5))