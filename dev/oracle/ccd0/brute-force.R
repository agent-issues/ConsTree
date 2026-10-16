# Brute-force validation of ConsTree::CCD0(): for n = 4..8 leaves, enumerate
# EVERY rooted binary tree, score it by sum log cc (plain-R clade counts), and
# check CCD0()'s tree, objective, normalizer, tree count and clade probabilities.
# usage: Rscript dev/oracle/ccd0/brute-force.R OUT.csv
source("dev/oracle/ccd0/referee.R")
library(ConsTree)
out <- commandArgs(TRUE)[1]

AllRootedTrees <- function(labels) {
  memo <- new.env()
  rec <- function(s) {
    key <- paste(s, collapse = "\r")
    if (!is.null(memo[[key]])) return(memo[[key]])
    if (length(s) == 1L) return(list(character(0)))
    rest <- s[-1L]
    res <- list()
    for (mask in 0:(2^length(rest) - 2)) {
      inA <- bitwAnd(mask, 2^(seq_along(rest) - 1)) > 0
      a <- c(s[1L], rest[inA]); b <- rest[!inA]
      ra <- rec(a); rb <- rec(b)
      for (ta in ra) for (tb in rb) res[[length(res) + 1L]] <- c(key, ta, tb)
    }
    memo[[key]] <- res
    res
  }
  lapply(rec(sort(labels)), function(cl) cl[grepl("\r", cl, fixed = TRUE)])
}

Sample <- function(type, n, k, labs) {
  switch(type,
    diffuse = lapply(seq_len(k), function(i) ape::rtree(n, tip.label = labs)),
    mixture = {
      base <- lapply(1:3, function(i) ape::rtree(n, tip.label = labs))
      w <- c(0.5, 0.3, 0.2)
      lapply(seq_len(k), function(i) {
        if (runif(1) < 0.8) base[[sample.int(3, 1, prob = w)]] else ape::rtree(n, tip.label = labs)
      })
    },
    few = lapply(seq_len(k), function(i) ape::rtree(n, tip.label = labs)),
    polytomy = lapply(seq_len(k), function(i) {
      tr <- ape::rtree(n, tip.label = labs)
      if (runif(1) < 0.3) {
        tr <- ape::di2multi(ape::compute.brlen(tr, function(m) sample(c(0, 1), m, TRUE, c(0.25, 0.75))))
        if (!TreeTools::TreeIsRooted(tr)) tr <- ape::rtree(n, tip.label = labs)
      }
      tr
    }))
}

rows <- list()
set.seed(20261009)
for (n in 4:8) {
  labs <- paste0("t", seq_len(n))
  cand <- AllRootedTrees(labs)
  stopifnot(length(cand) == prod(seq(1, 2 * n - 3, by = 2)))
  for (type in c("diffuse", "mixture", "few", "polytomy")) {
    for (rep in seq_len(if (n <= 6) 25 else 10)) {
      k <- switch(type, few = sample(2:4, 1), sample(c(10, 50, 200), 1))
      trees <- structure(Sample(type, n, k, sample(labs)), class = "multiPhylo")
      freq <- CladeFreq(trees)
      score <- vapply(cand, function(cl) {
        f <- unname(freq[cl]); if (anyNA(f)) -Inf else sum(log(f / k))
      }, double(1))
      best <- max(score)
      map <- tryCatch(CCD0(trees), error = function(e) e)
      if (inherits(map, "error")) {
        rows[[length(rows) + 1L]] <- data.frame(n, type, rep, k, nOpt = sum(score == best),
          status = if (!is.finite(best) && grepl("No fully resolved", conditionMessage(map))) "no-tree-ok" else "ERROR",
          objDiff = NA, inArgmax = NA, logPDiff = NA, nTreesOK = NA, probMaxDiff = NA)
        next
      }
      keys <- CladeKeys(map); keys <- keys[grepl("\r", keys, fixed = TRUE)]
      hit <- which(vapply(cand, function(cl) setequal(cl, keys), logical(1)))
      tol <- 1e-12 * (1 + abs(best))
      nOpt <- sum(score >= best - tol)
      fin <- is.finite(score)
      logZ <- log(sum(exp(score[fin] - best))) + best
      w <- exp(score - logZ)
      sub <- ape::subtrees(map)
      pBrute <- vapply(sub, function(s) {
        key <- paste(sort(s[["tip.label"]]), collapse = "\r")
        sum(w[fin & vapply(cand, function(cl) key %in% cl, logical(1))])
      }, double(1))
      rows[[length(rows) + 1L]] <- data.frame(n, type, rep, k, nOpt,
        status = "ok",
        objDiff = abs(attr(map, "logCredibility") - best),
        inArgmax = length(hit) == 1L && score[hit] >= best - tol,
        logPDiff = abs(attr(map, "logProbability") - (best - logZ)),
        nTreesOK = isTRUE(all.equal(exp(attr(map, "logNTrees")), sum(fin))),
        probMaxDiff = max(abs(attr(map, "probability") - pBrute)))
    }
    cat(n, type, "done\n")
  }
}
res <- do.call(rbind, rows)
write.csv(res, out, row.names = FALSE)
ok <- res$status == "ok"
cat("instances:", nrow(res), " computed:", sum(ok), " no-tree (correctly refused):",
    sum(res$status == "no-tree-ok"), " errors:", sum(res$status == "ERROR"), "\n")
cat("objective max |diff|:", max(res$objDiff[ok]), "\n")
cat("returned tree attains optimum:", sum(res$inArgmax[ok]), "/", sum(ok), "\n")
cat("unique optimum:", sum(res$nOpt[ok] == 1), "; tied optima:", sum(res$nOpt[ok] > 1), "\n")
cat("log P max |diff|:", max(res$logPDiff[ok]), "\n")
cat("tree count correct:", sum(res$nTreesOK[ok]), "/", sum(ok), "\n")
cat("clade probability max |diff|:", max(res$probMaxDiff[ok]), "\n")
print(aggregate(cbind(inArgmax, nOpt > 1) ~ n + type, data = res[ok, ], FUN = sum))
