# ===============================
# 1. Load data
# ===============================
GBM.vcf <- read.delim("data/GBM.vcf",
                      comment.char = "#",
                      header = TRUE,
                      stringsAsFactors = FALSE)

# ===============================
# 2. Basic summaries
# ===============================
cat("Total variants:", nrow(GBM.vcf), "\n")
cat("Variants per chromosome:\n")
print(table(GBM.vcf$CHROM))

# ===============================
# 3. Load gene annotation
# ===============================
ensembl71 <- read.delim("data/ensembl71.txt",
                        sep = "\t",
                        stringsAsFactors = FALSE,
                        header = TRUE)

# Merge gene info
GBM.vcf <- merge(ensembl71, GBM.vcf,
                 by.x = "ensg", by.y = "Gene",
                 all.x = FALSE, all.y = TRUE)

cat("Variants in genes:", nrow(GBM.vcf), "\n")

# ===============================
# 4. Filter somatic variants
# ===============================
GBM.somatic.vcf <- GBM.vcf[
  GBM.vcf$Normal.GT == "Homozygous.REF" &
  GBM.vcf$Normal.ALT.reads < 2, ]

cat("Somatic variants:", nrow(GBM.somatic.vcf), "\n")

# ===============================
# 5. Keep SNVs only
# ===============================
GBM.somatic.vcf <- GBM.somatic.vcf[
  nchar(GBM.somatic.vcf$REF) == 1 &
  nchar(GBM.somatic.vcf$ALT) == 1, ]

# ===============================
# 6. Fix chromosome names
# ===============================
if(length(grep("^chr", GBM.somatic.vcf$CHROM)) == 0) {
  GBM.somatic.vcf$CHROM <- paste0("chr", GBM.somatic.vcf$CHROM)
}

GBM.somatic.vcf <- GBM.somatic.vcf[
  GBM.somatic.vcf$CHROM %in% c(paste0("chr", 1:22), "chrX", "chrY"), ]

# ===============================
# 7. Load genome and extract sequence
# ===============================
library(BSgenome)
library(BSgenome.Hsapiens.UCSC.hg19)

GBM.somatic.vcf$upstream <- getSeq(Hsapiens,
                                  GBM.somatic.vcf$CHROM,
                                  start = GBM.somatic.vcf$POS - 1,
                                  width = 1,
                                  as.character = TRUE)

GBM.somatic.vcf$downstream <- getSeq(Hsapiens,
                                    GBM.somatic.vcf$CHROM,
                                    start = GBM.somatic.vcf$POS + 1,
                                    width = 1,
                                    as.character = TRUE)

# ===============================
# 8. Build trinucleotide context
# ===============================
GBM.somatic.vcf$context <- paste0(
  GBM.somatic.vcf$upstream,
  GBM.somatic.vcf$REF,
  GBM.somatic.vcf$downstream
)

GBM.somatic.vcf$mutation <- paste0(
  GBM.somatic.vcf$REF, "to", GBM.somatic.vcf$ALT
)

# ===============================
# 9. Normalize to 6 mutation types
# ===============================
# Convert A/G mutations to complementary (C/T framework)

complement <- function(base) {
  chartr("ACGT", "TGCA", base)
}

normalize_mut <- function(ref, alt, context) {
  if(ref %in% c("A", "G")) {
    ref <- complement(ref)
    alt <- complement(alt)
    context <- paste0(
      complement(substr(context,3,3)),
      ref,
      complement(substr(context,1,1))
    )
  }
  return(list(ref=ref, alt=alt, context=context))
}

norm <- mapply(normalize_mut,
               GBM.somatic.vcf$REF,
               GBM.somatic.vcf$ALT,
               GBM.somatic.vcf$context,
               SIMPLIFY = FALSE)

GBM.somatic.vcf$ref <- sapply(norm, `[[`, "ref")
GBM.somatic.vcf$alt <- sapply(norm, `[[`, "alt")
GBM.somatic.vcf$context <- sapply(norm, `[[`, "context")

GBM.somatic.vcf$mutation <- paste0(
  GBM.somatic.vcf$ref, "to", GBM.somatic.vcf$alt
)

# ===============================
# 10. Count frequencies
# ===============================
contexts <- unique(GBM.somatic.vcf$context)
mut.types <- c("CtoA","CtoG","CtoT","TtoA","TtoC","TtoG")

freq <- table(GBM.somatic.vcf$mutation, GBM.somatic.vcf$context)
freq <- prop.table(freq)

# ===============================
# 11. Plot signature
# ===============================
cols <- c("skyblue","black","red","grey","green","pink")

par(mfrow=c(1,6), mar=c(8,4,4,1))

for(i in 1:length(mut.types)) {
  mut <- mut.types[i]
  vals <- freq[mut, ]
  vals <- vals[order(names(vals))]
  
  barplot(vals,
          col = cols[i],
          border = "black",
          main = mut,
          las = 2,
          cex.names = 0.6,
          ylab = "frequency")
}