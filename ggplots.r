expression <- read.delim("data/expression.txt", stringsAsFactors = TRUE)

library(ggplot2)

gp <- ggplot(data = expression,
             mapping = aes(x = subtype, y = ESR1, fill = subtype))

# print(
#   gp + geom_boxplot()
# )

# print(gp + geom_violin(trim = FALSE))
# install.packages("ggbeeswarm")
library(ggbeeswarm)
print(gp + geom_quasirandom(size = 3))
print(gp + 
  geom_boxplot(outlier.shape = NA) +
  geom_quasirandom(size = 1.5, width = 0.25) +
  scale_fill_manual(values = c("red", "pink", "blue", "purple", "brown")))

patients <- read.delim("data/patients.txt", na.strings="", stringsAsFactors = FALSE)

patients_sub <- subset(patients, !is.na(subgroup))
gp <- ggplot(data = patients_sub,
  mapping = aes(x = location, fill = subgroup))

print(gp + geom_bar())

patients_sub$subgroup <- as.factor(patients_sub$subgroup)
levels(patients_sub$subgroup)

levels(patients_sub$subgroup) <- c("H3.3_K27M",
"H3.1_K27M","H3.3_G34R","BRAF","IDH1","wt")

glioma_colours <- c("H3.3_K27M" = "green", 
"H3.1_K27M" = "dark green", 
"H3.3_G34R" = "blue", 
"BRAF" = "gold", 
"IDH1" = "red",
"wt" = "grey")