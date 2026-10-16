# Build the rooted posterior samples used to compare Hipstr() / MrHipstr()
# with TreeAnnotator X and sumt.  Writes <outDir>/<name>.nex (BEAST-style NEXUS,
# Translate block, [&R]) for each data set, plus datasets.csv.
#
#   Rscript make-datasets.R <lib> <ebov.trees> <outDir>
#
# Real data: the 1 000-tree posterior sample of the 1 610-genome West African
# EBOV analysis of Dudas et al. (2017), Makona_1610_cds_ig.1K.trees, from
# https://github.com/ebov/space-time (Analyses/Phylogenetic) -- the 1 610-tip
# data set of Baele et al. (2025), Table 1.  Tip labels are recoded T1..T1610
# (first-tree order).  Smaller real samples prune the full posterior to random
# leaf subsets (ape::keep.tip, which preserves rooting and edge lengths).
#
# Simulated data: a coalescent tree (ape::rcoal) perturbed independently for
# each sample by Poisson(lambda) rooted NNI moves, edges chosen with
# probability inversely proportional to their length, so that short edges
# are poorly supported, as in real posteriors.

args <- commandArgs(TRUE)
.libPaths(c(args[[1]], .libPaths()))
ebovFile <- args[[2]]
outDir <- args[[3]]
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                   value = TRUE))), "referee.R"))

RootedNNI <- function(tr) {
  n <- length(tr[["tip.label"]])
  edge <- tr[["edge"]]
  internal <- setdiff(n + seq_len(tr[["Nnode"]]), n + 1L)
  len <- tr[["edge.length"]][match(internal, edge[, 2])]
  v <- internal[sample.int(length(internal), 1L,
                           prob = 1 / (len + 1e-3 * mean(len)))]
  p <- edge[edge[, 2] == v, 1]
  s <- edge[edge[, 1] == p & edge[, 2] != v, 2]
  kidsV <- edge[edge[, 1] == v, 2]
  x <- kidsV[[sample.int(2L, 1L)]]
  ex <- which(edge[, 1] == v & edge[, 2] == x)
  es <- which(edge[, 1] == p & edge[, 2] == s)
  edge[ex, 2] <- s
  edge[es, 2] <- x
  el <- tr[["edge.length"]]
  el[c(ex, es)] <- el[c(es, ex)]
  tr[["edge"]] <- edge
  tr[["edge.length"]] <- el
  attr(tr, "order") <- NULL
  tr
}

SimPosterior <- function(n, N, lambda, seed) {
  set.seed(seed)
  base <- ape::rcoal(n, tip.label = paste0("T", seq_len(n)))
  lapply(seq_len(N), function(i) {
    tr <- base
    for (j in seq_len(rpois(1L, lambda))) tr <- RootedNNI(tr)
    ape::reorder.phylo(tr, "cladewise")
  })
}

info <- list()
Save <- function(name, trees, source) {
  WriteBeastNexus(trees, file.path(outDir, paste0(name, ".nex")))
  info[[length(info) + 1L]] <<- data.frame(
    dataset = name, nTip = length(trees[[1]][["tip.label"]]),
    nTree = length(trees), source = source)
  cat(name, length(trees[[1]][["tip.label"]]), length(trees), "\n")
}

t0 <- proc.time()
ebov <- ReadTreesFile(ebovFile)
cat("Read EBOV:", length(ebov), "trees,", length(ebov[[1]][["tip.label"]]),
    "tips in", (proc.time() - t0)[["elapsed"]], "s\n")
orig <- ebov[[1]][["tip.label"]]
code <- setNames(paste0("T", seq_along(orig)), orig)
ebov <- lapply(ebov, function(tr) {
  tr[["tip.label"]] <- unname(code[tr[["tip.label"]]])
  tr
})
stopifnot(all(vapply(ebov, ape::is.binary, logical(1))),
          all(vapply(ebov, ape::is.rooted, logical(1))))
write.csv(data.frame(code = code, label = names(code)),
          file.path(outDir, "ebov_tip_codes.csv"), row.names = FALSE)
Save("ebov_n1610_k1000", ebov, "EBOV Dudas et al. 2017 (full)")
set.seed(2025)
for (n in c(400, 200, 100, 50, 20)) {
  keep <- sample(paste0("T", seq_along(orig)), n)
  sub <- lapply(ebov, ape::keep.tip, keep)
  Save(sprintf("ebov_n%d_k1000", n), sub, "EBOV pruned to random leaves")
  if (n == 100) {
    Save("ebov_n100_k100", sub[1:100], "EBOV pruned, first 100 trees")
  }
}
sims <- list(c(12, 200, 4), c(30, 5000, 3), c(50, 2000, 5), c(100, 5000, 8),
             c(200, 3000, 12), c(300, 1000, 15), c(400, 4000, 20))
for (i in seq_along(sims)) {
  s <- sims[[i]]
  Save(sprintf("sim_n%d_k%d", s[[1]], s[[2]]),
       SimPosterior(s[[1]], s[[2]], s[[3]], seed = 100 + i),
       sprintf("simulated: rcoal + Poisson(%g) rooted NNI", s[[3]]))
}
write.csv(do.call(rbind, info), file.path(outDir, "datasets.csv"),
          row.names = FALSE)
