# Real-data inputs from the Dudas et al. (2017) EBOV posterior (1,000 trees,
# 1,610 tips; github.com/ebov/space-time Analyses/Phylogenetic): random leaf
# subsets of every tree (the marginal posterior on those leaves).
# usage: Rscript ebov-prep.R EBOV.trees OUTDIR
args <- commandArgs(TRUE)
source("dev/oracle/ccd0/beast-io.R")
trees <- ape::.uncompressTipLabel(ReadBeastTrees(args[1]))
labels <- sort(trees[[1]][["tip.label"]])
cat("EBOV:", length(trees), "trees,", length(labels), "tips\n")
set.seed(1610)
for (n in c(50, 100, 200, 400, 800)) {
  tips <- sort(sample(labels, n))
  pruned <- lapply(trees, ape::keep.tip, tips)
  WriteBeastNexus(pruned, file.path(args[2], sprintf("EBOV_n%d_k1000.trees", n)), tips)
  WriteBeastNexus(pruned[seq(10, 1000, by = 10)],
                  file.path(args[2], sprintf("EBOV_n%d_k100.trees", n)), tips)
}
