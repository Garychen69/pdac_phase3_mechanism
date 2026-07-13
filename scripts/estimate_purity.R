# estimate_purity.R
# Compute real ESTIMATE (Yoshihara et al. 2013, Nat Commun) stromal/immune/purity
# scores for the three Phase 2 bulk cohorts, replacing the ad hoc 8-gene
# expression-mean purity proxy previously used in deconvolve_bulk_purity.py.
#
# Uses the `tidyestimate` CRAN package, which reimplements the published
# ESTIMATE algorithm (same gene sets: 141-gene stromal_signature and 141-gene
# immune_signature; same ssGSEA-style enrichment score computation) rather
# than the abandoned R-Forge `estimate` package, which is no longer installable.
#
# IMPORTANT platform caveat: the published purity-conversion formula
# (purity = cos(0.6049872018 + 0.0001467884 * estimate_score)) was calibrated
# on Affymetrix array data and is only valid for Affymetrix cohorts. Of the
# three Phase 2 bulk cohorts:
#   - GSE79668: RNA-seq (Illumina HiSeq)          -> formula NOT applied
#   - GSE71729: Agilent microarray                -> formula NOT applied
#   - GSE62165: Affymetrix HG-U219                 -> formula applied (valid)
# For the two non-Affymetrix cohorts we use the raw ESTIMATE score
# (stromal + immune enrichment) as a continuous stromal/immune-content
# covariate, which does not require the platform-specific conversion.

suppressMessages(library(tidyestimate))
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
expr_path <- args[1]
out_path <- args[2]
is_affy <- as.logical(args[3])

cat(sprintf("Reading expression matrix: %s\n", expr_path))
expr <- read.delim(gzfile(expr_path), header = TRUE, row.names = 1, check.names = FALSE)
cat(sprintf("  %d genes x %d samples\n", nrow(expr), ncol(expr)))

# Our expression matrices are already indexed by HGNC gene symbol, so no
# entrez/hgnc ID conversion is needed (filter_common_genes is for that step
# and, if used here, leaves a stray non-numeric "gene" id column behind that
# estimate_score mistakes for an extra sample). estimate_score() itself
# intersects rownames(df) against the stromal/immune gene sets, so passing
# the symbol-indexed matrix directly is sufficient and correct.
df <- data.frame(gene = rownames(expr), expr, check.names = FALSE)

scores <- estimate_score(df, is_affymetrix = is_affy)

if (!is_affy) {
  scores$purity <- NA_real_
}

write.table(scores, out_path, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf("Saved: %s (%d samples)\n", out_path, nrow(scores)))
