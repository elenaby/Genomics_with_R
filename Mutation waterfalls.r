somatic.variants <- read.delim("data/Lawrence.S2.txt", stringsAsFactors=T)
cancer_rates <- tapply(somatic.variants$logn_coding_mutations, somatic.variants$tumor_type, median)
cancer_rates <- cancer_rates[order(cancer_rates, decreasing=F)]
somatic.variants$tumor_type <- factor(somatic.variants$tumor_type, levels = names(cancer_rates))
somatic.variants <- somatic.variants[order(somatic.variants$tumor_type, somatic.variants$n_coding_mutations), ]
install.packages("GGally")
library(GGally)

ggplot(data = somatic.variants,
       mapping = aes(x = tumor_type,
                     y = log10(n_coding_mutations+1))) +
  #add alternating background colours
  geom_stripped_cols(odd = "#FFFFFF", even = "#CFE5DF") +
  #add out scatter points
  geom_point(position = position_dodge2(width = 0.9),
             size = 0.25) +
  #change Y axis labels
  scale_y_continuous(labels = c(0,10,100,1000,10000), expand = c(0,0)) +
  #add dotted lines
  geom_hline(yintercept = c(1,2,3,4), linetype = "dotted") +
  #change theme
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        panel.grid = element_blank()) +
  labs(y = "Coding mutations count") +
  #add a line for the median value
  stat_summary(fun = median,
               geom="crossbar",
               size = 0.25,
               width = 0.5,
               group = 1,
               show.legend = FALSE,
               color = "#000000")