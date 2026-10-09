# Exhaustive check of Hipstr() / MrHipstr() on small inputs: enumerate every
# rooted binary tree on n <= 8 leaves (135 135 trees at n = 8), score each
# under the HIPSTR objective with an independent bitmask clade table, and
# confirm that the returned tree attains the optimum.
#
# Run (Hamilton):  Rscript dev/oracle/hipstr/brute-force.R <lib> <outfile>

args <- commandArgs(TRUE)
if (length(args) >= 1) .libPaths(c(args[[1]], .libPaths()))
outFile <- if (length(args) >= 2) args[[2]] else "brute-force.csv"
library(ConsTree)

# --- every rooted binary tree on leaf set S (bitmask), memoised -------------
# Returns list(clades = matrix [nTree x (|S|-1)] of clade masks,
#              keys   = matrix of pair keys C * 2^n + A, A the child holding
#                       the lowest leaf of C)
.memo <- new.env()
Enumerate <- function(S, n) {
  k <- paste(n, S)
  if (!is.null(.memo[[k]])) return(.memo[[k]])
  bits <- which(bitwAnd(S, bitwShiftL(1L, 0:(n - 1L))) != 0L) - 1L
  if (length(bits) == 1L) {
    res <- list(clades = matrix(integer(0), 1, 0), keys = matrix(numeric(0), 1, 0))
  } else {
    low <- bitwShiftL(1L, bits[[1]])
    rest <- bits[-1]
    cl <- list()
    ky <- list()
    for (m in 0:(2^length(rest) - 2)) {          # proper subsets of rest
      A <- low
      for (j in seq_along(rest)) {
        if (bitwAnd(m, bitwShiftL(1L, j - 1L))) A <- bitwOr(A, bitwShiftL(1L, rest[[j]]))
      }
      B <- bitwXor(S, A)
      ta <- Enumerate(A, n)
      tb <- Enumerate(B, n)
      ia <- rep(seq_len(nrow(ta$clades)), each = nrow(tb$clades))
      ib <- rep(seq_len(nrow(tb$clades)), times = nrow(ta$clades))
      cl[[length(cl) + 1L]] <- cbind(ta$clades[ia, , drop = FALSE],
                                     tb$clades[ib, , drop = FALSE], S)
      ky[[length(ky) + 1L]] <- cbind(ta$keys[ia, , drop = FALSE],
                                     tb$keys[ib, , drop = FALSE],
                                     S * 2^n + A)
    }
    res <- list(clades = do.call(rbind, cl), keys = do.call(rbind, ky))
  }
  .memo[[k]] <- res
  res
}

# Clade masks of a phylo, as (clade, child-with-lowest-leaf) per internal node
TreeMasks <- function(tr, labels) {
  n <- length(labels)
  tr <- TreeTools::RenumberTips(tr, labels)
  po <- ape::reorder.phylo(tr, "postorder")[["edge"]]
  mask <- numeric(n + tr[["Nnode"]])
  mask[seq_len(n)] <- 2^(seq_len(n) - 1)
  lowChild <- numeric(n + tr[["Nnode"]])
  low <- c(seq_len(n), rep(Inf, tr[["Nnode"]]))
  for (r in seq_len(nrow(po))) {
    p <- po[r, 1]; ch <- po[r, 2]
    mask[p] <- mask[p] + mask[ch]
    if (low[ch] < low[p]) { low[p] <- low[ch]; lowChild[p] <- mask[ch] }
  }
  internal <- n + seq_len(tr[["Nnode"]])
  cbind(clade = mask[internal], key = mask[internal] * 2^n + lowChild[internal])
}

Score <- function(trees, labels) {
  n <- length(labels)
  N <- length(trees)
  count <- numeric(2^n)
  pairs <- numeric(0)
  for (tr in trees) {
    m <- TreeMasks(tr, labels)
    count[m[, "clade"]] <- count[m[, "clade"]] + 1
    pairs <- c(pairs, m[, "key"])
  }
  list(count = count, pairs = unique(pairs), N = N)
}

ScoreTrees <- function(en, tab) {
  cnt <- matrix(tab$count[en$clades], nrow(en$clades))
  ok <- matrix(en$keys %in% tab$pairs, nrow(en$keys))
  valid <- rowSums(!ok) == 0
  list(valid = valid,
       logCC = rowSums(log(cnt / tab$N)),
       nMaj = rowSums(2 * cnt > tab$N))
}

RandomSample <- function(n, N, nBase, seed, mode) {
  set.seed(seed)
  labs <- paste0("t", seq_len(n))
  if (mode == "mixture") {
    base <- lapply(seq_len(nBase), function(i) TreeTools::RandomTree(labs, root = TRUE))
    return(base[sample.int(nBase, N, replace = TRUE)])
  }
  # 'uniform': independent uniformly random rooted trees (very diffuse)
  lapply(seq_len(N), function(i) TreeTools::RandomTree(labs, root = TRUE))
}

results <- list()
cfgs <- expand.grid(n = 3:8, N = c(4, 10, 25, 100), mode = c("mixture", "uniform"),
                    stringsAsFactors = FALSE)
for (i in seq_len(nrow(cfgs))) {
  n <- cfgs$n[[i]]; N <- cfgs$N[[i]]; mode <- cfgs$mode[[i]]
  nRep <- if (n <= 6) 25 else if (n == 7) 40 else 25
  labs <- paste0("t", seq_len(n))
  en <- Enumerate(2^n - 1, n)
  for (rep in seq_len(nRep)) {
    seed <- 1e5 * n + 1e3 * N + rep + if (mode == "uniform") 5e7 else 0
    trees <- RandomSample(n, N, nBase = max(2, round(N / 3)), seed, mode)
    tab <- Score(trees, labs)
    sc <- ScoreTrees(en, tab)
    bestHip <- max(sc$logCC[sc$valid])
    nOptHip <- sum(sc$valid & abs(sc$logCC - bestHip) < 1e-9)
    bestMaj <- max(sc$nMaj[sc$valid])
    bestMr <- max(sc$logCC[sc$valid & sc$nMaj == bestMaj])
    nOptMr <- sum(sc$valid & sc$nMaj == bestMaj & abs(sc$logCC - bestMr) < 1e-9)

    for (fn in c("Hipstr", "MrHipstr")) {
      res <- get(fn)(structure(trees, class = "multiPhylo"))
      m <- TreeMasks(res, labs)
      cnt <- tab$count[m[, "clade"]]
      inSpace <- all(m[, "key"] %in% tab$pairs)
      lc <- sum(log(cnt / N))
      nm <- sum(2 * cnt > N)
      ok <- if (fn == "Hipstr") {
        inSpace && abs(lc - bestHip) < 1e-9
      } else {
        inSpace && nm == bestMaj && abs(lc - bestMr) < 1e-9
      }
      results[[length(results) + 1L]] <- data.frame(
        method = fn, n = n, N = N, mode = mode, rep = rep, seed = seed,
        nCandidates = nrow(en$clades), nInSpace = sum(sc$valid),
        optimum = if (fn == "Hipstr") bestHip else bestMr,
        optMaj = if (fn == "Hipstr") NA else bestMaj,
        nOptimal = if (fn == "Hipstr") nOptHip else nOptMr,
        returned = lc, returnedMaj = nm, inSpace = inSpace,
        attrMatches = abs(attr(res, "logCredibility") - lc) < 1e-12,
        optimal = ok)
    }
  }
}
out <- do.call(rbind, results)
write.csv(out, outFile, row.names = FALSE)
cat("datasets x methods:", nrow(out), " optimal:", sum(out$optimal),
    " in space:", sum(out$inSpace), " attr ok:", sum(out$attrMatches),
    " with ties (>1 optimal tree):", sum(out$nOptimal > 1), "\n")
print(aggregate(optimal ~ method + n, out, function(x) paste0(sum(x), "/", length(x))))
