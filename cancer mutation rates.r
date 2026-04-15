mutation.rates <- read.delim("data/Vogelstein.mutation.rates.txt", stringsAsFactors=T)
barplot(mutation.rates$Average, names.arg= mutation.rates$Tumor.Type)
mutation.rates <- mutation.rates[order(mutation.rates$Average),]
barplot(mutation.rates$Average, names.arg= mutation.rates$Tumor.Type, las=2)
barplot(mutation.rates$Average, names.arg= mutation.rates$Tumor.Type, las=2, col=rainbow(nrow(mutation.rates)), cex.names=0.75, log="y")