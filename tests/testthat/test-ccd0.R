# Rooted clades of a tree as sorted label strings (size 2..n-1, plus the root)
ccd0Clades <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  edge <- tree[["edge"]]
  nTip <- TreeTools::NTip(tree)
  tl <- tree[["tip.label"]]
  desc <- vector("list", nTip + tree[["Nnode"]])
  for (i in seq_len(nTip)) desc[[i]] <- tl[i]
  for (e in rev(seq_len(nrow(edge)))) {
    desc[[edge[e, 1]]] <- c(desc[[edge[e, 1]]], desc[[edge[e, 2]]])
  }
  keys <- vapply(desc[-seq_len(nTip)], function(x) paste(sort(x), collapse = ","),
                 character(1))
  sort(unique(keys))
}

# Independent brute-force referee: every rooted binary tree on `labels`, as a
# set of clade keys (memoized recursion on the leaf set).
allRootedTrees <- function(labels) {
  memo <- new.env()
  rec <- function(s) {
    key <- paste(s, collapse = ",")
    if (!is.null(memo[[key]])) return(memo[[key]])
    if (length(s) == 1L) return(list(character(0)))
    rest <- s[-1L]
    out <- list()
    # Child A contains s[1]; B is the non-empty remainder
    for (mask in 0:(2^length(rest) - 2)) {
      inA <- bitwAnd(mask, 2^(seq_along(rest) - 1)) > 0
      a <- c(s[1L], rest[inA])
      b <- rest[!inA]
      for (ta in rec(a)) for (tb in rec(b)) {
        out[[length(out) + 1L]] <- c(key, ta, tb)
      }
    }
    memo[[key]] <- out
    out
  }
  lapply(rec(sort(labels)), function(cl) {
    cl[lengths(strsplit(cl, ",")) > 1L]
  })
}

bruteCCD0 <- function(trees) {
  labels <- sort(TreeTools::TipLabels(trees[[1]]))
  k <- length(trees)
  freq <- table(unlist(lapply(trees, ccd0Clades)))
  cand <- allRootedTrees(labels)
  score <- vapply(cand, function(cl) {
    f <- as.numeric(freq[cl])
    if (anyNA(f)) -Inf else sum(log(f / k))
  }, double(1))
  list(cand = cand, score = score)
}

test_that("CCD0 reproduces the worked example of Berling et al. (2025, fig. 1)", {
  trees <- ape::read.tree(text = c(
    rep("(((A, B), C), (D, E));", 3),
    rep("(((A, (B, C)), D), E);", 2),
    rep("(((A, (B, C)), E), D);", 2)
  ))
  map <- CCD0(trees)
  expect_s3_class(map, "phylo")
  expect_true(ape::is.rooted(map))
  expect_true(ape::is.binary(map))
  # The MAP tree was never sampled
  expect_equal(ccd0Clades(map),
               ccd0Clades(ape::read.tree(text = "((A, (B, C)), (D, E));")))
  expect_equal(exp(attr(map, "logProbability")), 12 / 49)
  expect_equal(exp(attr(map, "logCredibility")), 12 / 49)
  expect_equal(exp(attr(map, "logNTrees")), 6)
  expect_equal(attr(map, "nClades"), 12L)  # 5 leaves, root, 6 clades
  expect_equal(attr(map, "nSplits"), 10)
  # Root, then clades in preorder
  expect_equal(sort(attr(map, "credibility")), sort(c(1, 1, 4 / 7, 3 / 7)))
  expect_equal(as.numeric(map[["node.label"]]), attr(map, "credibility"))
  expect_equal(attr(map, "probability")[1], 1)
  expect_identical(CCD0MAP(trees), map)
})

test_that("CCD0 includes unobserved splits and breaks ties towards balance", {
  # Fig. 2 of Berling et al.: clades AB and CD are observed, split AB|CD is not.
  trees <- ape::read.tree(text = c("(((A, B), C), D);", "(((C, D), A), B);"))
  map <- CCD0(trees)
  # Three splits of the root tie at 1/4; the balanced one wins
  expect_equal(ccd0Clades(map),
               ccd0Clades(ape::read.tree(text = "((A, B), (C, D));")))
  expect_equal(exp(attr(map, "logProbability")), 1 / 3)
  expect_equal(exp(attr(map, "logNTrees")), 3)
})

test_that("CCD0 ties of equal balance go to the lexicographically first child", {
  # Root splits {A}|{BCD} and {D}|{ABC}: both size-1 smaller children; A first.
  trees <- ape::read.tree(text = c("(A, ((B, C), D));", "(((B, C), A), D);"))
  map <- CCD0(trees)
  expect_equal(ccd0Clades(map),
               ccd0Clades(ape::read.tree(text = "(A, ((B, C), D));")))
  # "First" is alphabetical, whatever the tip numbering or tree order
  rev1 <- TreeTools::RenumberTips(trees[[1]], rev(TreeTools::TipLabels(trees[[1]])))
  expect_equal(ccd0Clades(CCD0(c(rev1, trees[[2]]))), ccd0Clades(map))
  expect_equal(ccd0Clades(CCD0(c(trees[[2]], rev1))), ccd0Clades(map))
  # ... compared character by character (radix order: "B" < "a")
  relab <- lapply(trees, function(tr) {
    tr[["tip.label"]] <- c(A = "a", B = "c", C = "d", D = "B")[tr[["tip.label"]]]
    tr
  })
  expect_equal(ccd0Clades(CCD0(relab)),
               ccd0Clades(ape::read.tree(text = "(((c, d), a), B);")))
})

test_that("CCD0 attains the brute-force optimum", {
  set.seed(2817)
  labs <- letters[1:6]
  base <- lapply(1:3, function(i) ape::rtree(6, tip.label = labs))
  for (rep in 1:4) {
    trees <- structure(lapply(seq_len(25), function(i) {
      if (runif(1) < 0.6) base[[sample.int(3, 1)]] else ape::rtree(6, tip.label = sample(labs))
    }), class = "multiPhylo")
    map <- CCD0(trees)
    bf <- bruteCCD0(trees)
    best <- max(bf[["score"]])
    expect_equal(attr(map, "logCredibility"), best, tolerance = 1e-12)
    mapKeys <- ccd0Clades(map)
    hit <- vapply(bf[["cand"]], function(cl) setequal(cl, mapKeys), logical(1))
    expect_equal(sum(hit), 1L)
    expect_equal(bf[["score"]][hit], best, tolerance = 1e-12)
    finite <- is.finite(bf[["score"]])
    logZ <- log(sum(exp(bf[["score"]][finite])))
    expect_equal(attr(map, "logProbability"), best - logZ, tolerance = 1e-10)
    expect_equal(exp(attr(map, "logNTrees")), sum(finite))
    # CCD0 clade probabilities = posterior mass of trees containing the clade
    w <- exp(bf[["score"]] - logZ)
    sub <- ape::subtrees(map)
    for (i in seq_along(sub)) {
      key <- paste(sort(sub[[i]][["tip.label"]]), collapse = ",")
      pBrute <- sum(w[finite & vapply(bf[["cand"]], function(cl) key %in% cl,
                                       logical(1))])
      expect_equal(attr(map, "probability")[i], pBrute, tolerance = 1e-10)
    }
  }
})

test_that("CCD0 tree is at least as credible as every sampled tree", {
  set.seed(5021)
  labs <- paste0("t", 1:30)
  base <- ape::rtree(30, tip.label = labs)
  trees <- structure(lapply(1:60, function(i) {
    tr <- base
    # perturb by swapping two leaf labels
    sw <- sample(30, 2)
    tr[["tip.label"]][sw] <- tr[["tip.label"]][rev(sw)]
    tr
  }), class = "multiPhylo")
  map <- CCD0(trees)
  freq <- table(unlist(lapply(trees, ccd0Clades)))
  sampled <- vapply(trees, function(tr) sum(log(freq[ccd0Clades(tr)] / 60)),
                    double(1))
  expect_gte(attr(map, "logCredibility"), max(sampled) - 1e-12)
  mapScore <- sum(log(as.numeric(freq[ccd0Clades(map)]) / 60))
  expect_equal(mapScore, attr(map, "logCredibility"), tolerance = 1e-12)
  expect_lte(attr(map, "logProbability"), 0)
  expect_true(all(attr(map, "probability") <= 1 + 1e-12))
  expect_true(all(attr(map, "probability") > 0))
})

test_that("CCD0 is idempotent and invariant to tip numbering", {
  tree <- ape::rtree(12)
  map <- CCD0(c(tree, tree, tree))
  expect_equal(ccd0Clades(map), ccd0Clades(tree))
  expect_equal(attr(map, "credibility"), rep(1, 11))
  expect_equal(attr(map, "logProbability"), 0)
  expect_equal(attr(map, "probability"), rep(1, 11))
  single <- CCD0(tree)
  expect_equal(ccd0Clades(single), ccd0Clades(tree))

  set.seed(77)
  trees <- structure(lapply(1:20, function(i) ape::rtree(9, tip.label = letters[1:9])),
                     class = "multiPhylo")
  ref <- CCD0(trees)
  shuffled <- structure(lapply(seq_along(trees), function(i) {
    if (i == 1) trees[[1]] else
      TreeTools::RenumberTips(trees[[i]], sample(letters[1:9]))
  }), class = "multiPhylo")
  expect_equal(ccd0Clades(CCD0(shuffled)), ccd0Clades(ref))
  expect_equal(attr(CCD0(shuffled), "logProbability"), attr(ref, "logProbability"))
  # Compressed tip labels
  compressed <- ape::.compressTipLabel(trees)
  expect_equal(ccd0Clades(CCD0(compressed)), ccd0Clades(ref))
  # Plain list with NULL entries
  expect_equal(ccd0Clades(CCD0(c(list(NULL), unclass(trees)))), ccd0Clades(ref))
  noProb <- CCD0(trees, cladeProbability = FALSE)
  expect_true(all(is.na(attr(noProb, "probability"))))
  expect_equal(ccd0Clades(noProb), ccd0Clades(ref))
  # Output is a valid phylo object
  expect_equal(ape::write.tree(ape::read.tree(text = ape::write.tree(ref))),
               ape::write.tree(ref))
  expect_equal(TreeTools::NTip(ref), 9L)
  expect_equal(nrow(ref[["edge"]]), 16L)
})

test_that("CCD0 handles small and polytomous inputs", {
  two <- ape::read.tree(text = c("(A, B);", "(B, A);"))
  map2 <- CCD0(two)
  expect_equal(TreeTools::NTip(map2), 2L)
  expect_equal(map2[["Nnode"]], 1L)
  expect_equal(attr(map2, "logProbability"), 0)

  poly <- ape::read.tree(text = c("((A, B, C), D);", "((A, B, C), D);",
                                  "(((A, B), C), D);"))
  map <- CCD0(poly)
  expect_equal(ccd0Clades(map),
               ccd0Clades(ape::read.tree(text = "(((A, B), C), D);")))
  expect_equal(attr(map, "credibility"), c(1, 1, 1 / 3))

  star <- ape::read.tree(text = c("((A, B, C), D);", "((A, B, C), D);"))
  expect_error(CCD0(star), "No fully resolved tree")

  # Singleton nodes do not double-count a clade
  single <- ape::read.tree(text = "((((A, B)), C), D);")
  mapS <- CCD0(c(single, ape::read.tree(text = "(((A, B), C), D);")))
  expect_equal(attr(mapS, "credibility"), c(1, 1, 1))
})

test_that("CCD0 rejects invalid input", {
  expect_error(CCD0(1), "must be a list")
  expect_error(CCD0(list()), "no trees")
  expect_error(CCD0(list(ape::rtree(4), 1)), "class `phylo`")
  expect_error(CCD0(c(ape::rtree(4), ape::rtree(5))), "same tip labels")
  dup <- ape::read.tree(text = "((A, A), B);")
  expect_error(CCD0(c(dup, dup)), "unique")
  expect_error(CCD0(ape::read.tree(text = c("(A, B, C);", "((A, B), C);"))),
               "unrooted")
  many <- ape::read.tree(text = rep("(A, B, C);", 7))
  expect_error(CCD0(many), "1, 2, 3, 4, 5, \\.\\.\\.")
  one <- list(structure(list(edge = matrix(c(2L, 1L), 1), tip.label = "A",
                             Nnode = 1L), class = "phylo"))
  expect_error(CCD0(one), "at least two leaves")
})
