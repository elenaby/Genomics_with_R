CNA <- read.delim("data/CNA.txt")
plot(CNA$Amplification)
CNA <- CNA[order(CNA$chrom, CNA$start), ]
head(CNA)
plot(CNA$Amplification, type="h", col="red",
     main="Amplifications", ylab="frequency")

plot(-CNA$Deletion, type="h", col="green",
     main="Deletions", ylab="frequency")

CNA[CNA$Amplification > 0.5,]

CNA_amps <- CNA[CNA$Amplification > 0.3,]
write.table(CNA_amps, file="results/Common_Amps.xls", row.names=F, sep="\t")
CNA_dels <- CNA[CNA$Deletion > 0.2,]
write.table(CNA_dels, file="results/Common_Dels.xls", row.names=F, sep="\t")
CNA_common <- CNA[CNA$Amplification > 0.3|CNA$Deletion > 0.2,]
write.table(CNA_common, file="results/Common_CNAs.xls", row.names=F, sep="\t")