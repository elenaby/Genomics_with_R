BRCA_mutations <- read.delim(paste0("data/TCGA/BRCA.mutations.txt"))
BRCA_mutations$patientID <- gsub("-01A\\S+", "", BRCA_mutations$Tumor_Sample_Barcode)
BRCA_clin$TP53 <- "WT"
BRCA_clin$TP53[!BRCA_clin$barcode %in% BRCA_mutations$patientID] <- NA
BRCA_clin$TP53[
  BRCA_clin$barcode %in% 
  BRCA_mutations$patientID[BRCA_mutations$Hugo_Symbol == "TP53"]
] <- "MUT"
sfit_TP53 <- survfit(Surv(OS, status)~TP53, data=BRCA_clin)
coxph(Surv(OS, status)~TP53, data= BRCA_clin)
plot(sfit_TP53, col=c("black", "blue"), las=1, main="Breast cancer survival - TP53", xlab="Days")
legend("bottomleft", c( "TP53_MUT", "TP53_WT"), fill=c("black", "blue"), bty="n")
ggsurvplot(sfit_TP53, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
legend.labs=c("MUT", "WT"), legend.title="TP53",  
palette=c("black", "blue"), xlab="Days",
title="Kaplan-Meier Curve for Breast Cancer Survival", 
risk.table.height=.25)