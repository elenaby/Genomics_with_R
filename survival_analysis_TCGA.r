BRCA_clin <- read.delim("data/tcga_brca_clinical.txt",
 stringsAsFactor=T, na.strings=c("", "NA", " "), row.names=1)

sfit_ER <- survfit(Surv(OS, status)~ER, data=BRCA_clin)
coxph(Surv(OS, status)~ER, data= BRCA_clin)
plot(sfit_ER, col=c("black", "blue"), las=1, main="Breast cancer survival - ER", xlab="Days")
legend("bottomleft", c( "ER_negative", "ER_positive"), fill=c("black", "blue"), bty="n")
ggsurvplot(sfit_ER, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
  legend.labs=c("negative", "positive"), legend.title="ER",  
  palette=c("black", "blue"), xlab="Days",
  title="Kaplan-Meier Curve for Breast Cancer Survival", 
  risk.table.height=.25)