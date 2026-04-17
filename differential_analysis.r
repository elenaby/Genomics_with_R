library(Biobase)
NKI295 <- readExpressionSet(exprsFile="data/NKI295.exprs.txt", sep="\t", header=T, stringsAsFactors=F, row.names=1)
fData(NKI295) <- read.delim("data/NKI295.fdata.txt", sep="\t", header=T, row.names=1, stringsAsFactors=F)
pData(NKI295) <- read.table("data/NKI295.pdata.txt", sep="\t", header=T, row.names=1, stringsAsFactors=T) 
sampleNames(NKI295)[1:10]
exprs(NKI295)[c("ESR1", "ERBB2", "GATA3"),c("NKI295_12", "NKI295_123", "NKI295_268")]
subtype.colours <- c("red", "deeppink","dark blue", "light blue","orange")
boxplot(exprs(NKI295)[which(fData(NKI295)$symbol == "ESR1"),]~NKI295$subtype, col=subtype.colours)
NKI295 <- NKI295[,order(exprs(NKI295)["ESR1",])]
barplot(exprs(NKI295)["ESR1",], cex.names=0.1, las=2,col=subtype.colours[NKI295$subtype], names.arg="", border=NA)