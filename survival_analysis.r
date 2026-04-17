library(survival)
browseVignettes("survival")
sfit <- survfit(Surv(time, status)~sex, data=lung)
plot(sfit, col=c("blue", "pink"), main="Lung cancer")
survreg(Surv(time, status)~sex, data=lung)
coxph(Surv(time, status)~sex, data=lung)
coxph(Surv(time, status)~sex, data=lung)
library(survminer)
ggsurvplot(sfit, conf.int=TRUE, pval=TRUE, risk.table=TRUE, 
  legend.labs=c("Male", "Female"), legend.title="Sex",  
  palette=c("dodgerblue2", "orchid2"), 
  title="Kaplan-Meier Curve for Lung Cancer Survival", 
  risk.table.height=.15, xlab="Days")