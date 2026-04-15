# Create results folder if it doesn't exist
dir.create("results", showWarnings = FALSE)

# Load data
somatic.variants <- read.delim("data/Lawrence.S2.txt", stringsAsFactors = TRUE)

# ---- Step 1: Compute median mutation rates per tumor type ----
cancer_rates <- tapply(
  somatic.variants$logn_coding_mutations,
  somatic.variants$tumor_type,
  median,
  na.rm = TRUE
)

# Sort tumor types by median mutation rate (ascending)
cancer_rates <- cancer_rates[order(cancer_rates, decreasing = FALSE)]

# ---- Step 2: Reorder factor levels ----
somatic.variants$tumor_type <- factor(
  somatic.variants$tumor_type,
  levels = names(cancer_rates)
)

# ---- Step 3: Open PDF device ----
pdf("results/mutation_boxplot.pdf", width = 14, height = 8)

# Adjust margins for long labels
par(mar = c(12, 4, 4, 2))

# ---- Step 4: Boxplot ----
boxplot(logn_coding_mutations ~ tumor_type,
        data = somatic.variants,
        col = rainbow(length(cancer_rates)),
        pch = 16,
        cex = 0.5,
        las = 2,                # rotate labels
        xlab = "",
        main = "Cancer Mutation Rates",
        ylab = "Log Number of Somatic Mutations")

# ---- Step 5: Add stripchart (jittered points) ----
stripchart(logn_coding_mutations ~ tumor_type,
           data = somatic.variants,
           vertical = TRUE,
           method = "jitter",
           pch = 16,
           col = rgb(0, 0, 0, 0.3),  # semi-transparent points
           add = TRUE)

# ---- Step 6: Close PDF ----
dev.off()