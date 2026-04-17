library(limma)
NKI295$Basal <- ifelse(NKI295$subtype == "Basal", yes="Basal", no="Other")
NKI295$Basal <- factor(NKI295$Basal)
limma.parameters <- cbind(intercept=rep(1, ncol(NKI295)), NKI295$Basal)
limma.fit <- lmFit(exprs(NKI295), limma.parameters)
limma.fit <- eBayes(limma.fit)
limma.table <- topTable(limma.fit, coef=1, number=nrow(NKI295), sort.by="none")
limma.table <- cbind(limma.table, fData(NKI295))
Basal.significant <- limma.table[which(limma.table$adj.P.Val < 1e-4),]
Basal.significant <- Basal.significant[abs(Basal.significant$logFC) > 0.3,]
Basal.significant <- Basal.significant[order(Basal.significant$logFC, decreasing=T),]