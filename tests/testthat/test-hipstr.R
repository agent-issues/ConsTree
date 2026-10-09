# Independent referee for HIPSTR / MrHIPSTR, written without reference to the
# C++.  Clades are keyed by their sorted tip labels; a tree is in the HIPSTR
# search space iff each internal node splits its clade into a pair of child
# clades that occur as siblings in some input tree.

hipTriples <- function(tr) {
  n <- length(tr[["tip.label"]])
  edge <- tr[["edge"]]
  desc <- vector("list", n + tr[["Nnode"]])
  for (i in seq_len(n)) desc[[i]] <- tr[["tip.label"]][[i]]
  po <- ape::reorder.phylo(tr, "postorder")[["edge"]]
  for (r in seq_len(nrow(po))) {
    desc[[po[r, 1]]] <- c(desc[[po[r, 1]]], desc[[po[r, 2]]])
  }
  key <- function(v) paste(sort(desc[[v]]), collapse = ",")
  kids <- split(edge[, 2], edge[, 1])
  internal <- n + seq_len(tr[["Nnode"]])
  t(vapply(internal, function(v) {
    k <- kids[[as.character(v)]]
    c(key(v), sort(c(key(k[[1]]), key(k[[2]]))))
  }, character(3)))
}

hipTable <- function(trees) {
  trip <- lapply(trees, hipTriples)
  clades <- unlist(lapply(trip, function(x) x[, 1]))
  list(count = table(clades), nTree = length(trees),
       pairs = unique(unlist(lapply(trip, function(x)
         paste(x[, 1], x[, 2], x[, 3], sep = "|")))))
}

# c(nMajority, logCC), or NULL if the tree lies outside the search space
hipScore <- function(trip, tab) {
  if (!all(paste(trip[, 1], trip[, 2], trip[, 3], sep = "|") %in% tab$pairs)) {
    return(NULL)
  }
  cnt <- as.numeric(tab$count[trip[, 1]])
  c(sum(2 * cnt > tab$nTree), sum(log(cnt / tab$nTree)))
}

# Every rooted binary tree on `labs`, as a matrix of (clade, child, child) keys
allTriples <- function(labs) {
  labs <- sort(labs)
  if (length(labs) == 1L) return(list(matrix(character(0), 0, 3)))
  key <- function(x) paste(sort(x), collapse = ",")
  rest <- labs[-1]
  out <- list()
  for (r in 0:(length(rest) - 1L)) {
    combs <- if (r == 0L) list(character(0)) else
      utils::combn(rest, r, simplify = FALSE)
    for (cmb in combs) {
      a <- c(labs[[1]], cmb)
      b <- setdiff(rest, cmb)
      here <- c(key(labs), sort(c(key(a), key(b))))
      for (ta in allTriples(a)) for (tb in allTriples(b)) {
        out[[length(out) + 1L]] <- rbind(ta, tb, here)
      }
    }
  }
  out
}

bruteOptimum <- function(trees, majority) {
  tab <- hipTable(trees)
  scores <- lapply(allTriples(TreeTools::TipLabels(trees[[1]])), hipScore, tab)
  scores <- do.call(rbind, scores[!vapply(scores, is.null, logical(1))])
  if (majority) {
    m <- max(scores[, 1])
    c(m, max(scores[scores[, 1] == m, 2]))
  } else {
    c(NA, max(scores[, 2]))
  }
}

cladeKeys <- function(tr) sort(hipTriples(tr)[, 1])
nwkClades <- function(txt) cladeKeys(ape::read.tree(text = txt))

randomSample <- function(nTip, nTree, nBase, seed) {
  set.seed(seed)
  labs <- paste0("t", seq_len(nTip))
  base <- lapply(seq_len(nBase), function(i)
    TreeTools::RandomTree(labs, root = TRUE))
  structure(base[sample.int(nBase, nTree, replace = TRUE)],
            class = "multiPhylo")
}

test_that("Hipstr() reproduces a hand-checked example", {
  # 3 x ((a,d),(b,c)); 2 x ((a,b),(c,d)); 2 x (((a,b),c),d)
  trees <- ape::read.tree(text = c(
    rep("((a, d), (b, c));", 3),
    rep("((a, b), (c, d));", 2),
    rep("(((a, b), c), d);", 2)
  ))
  hip <- Hipstr(trees)
  expect_s3_class(hip, "phylo")
  expect_true(ape::is.rooted(hip))
  expect_true(ape::is.binary(hip))
  expect_equal(cladeKeys(hip), nwkClades("((a, d), (b, c));"))
  expect_equal(attr(hip, "logCredibility"), 2 * log(3 / 7))
  expect_equal(sort(attr(hip, "credibility")), c(3, 3, 7) / 7)
  expect_equal(hip[["node.label"]], as.character(attr(hip, "credibility")))
  expect_equal(attr(hip, "credibility")[[1]], 1)

  # The majority clade (a, b) [4 / 7] is absent from HIPSTR, forced by
  # MrHIPSTR.  ((a,b),(c,d)) and (((a,b),c),d) tie, at log(4/7) + log(2/7);
  # the larger summed child count (2 + 7 > 4 + 2) chooses the second.
  expect_false("a,b" %in% cladeKeys(hip))
  mr <- MrHipstr(trees)
  expect_equal(cladeKeys(mr), nwkClades("(((a, b), c), d);"))
  expect_equal(attr(mr, "logCredibility"), log(4 / 7) + log(2 / 7))
})

test_that("MrHipstr() rewards only clades in more than half the trees", {
  # (a, c) occurs in exactly 4 / 8 trees; (b, e) in 5 / 8.  Treating 1/2 as
  # 'majority' would give (((a, c), (b, e)), d) instead.
  trees <- ape::read.tree(text = c(
    rep("((a, c), (b, (d, e)));", 2),
    rep("((b, e), (c, (d, a)));", 3),
    "(((b, d), (a, e)), c);",
    rep("(d, ((a, c), (e, b)));", 2)
  ))
  expected <- nwkClades("(((a, d), c), (b, e));")
  expect_equal(cladeKeys(MrHipstr(trees)), expected)
  expect_equal(cladeKeys(Hipstr(trees)), expected)
  expect_false("a,c" %in% cladeKeys(MrHipstr(trees)))
})

test_that("ties are broken by the documented rule, independent of order", {
  # Each three-leaf resolution once: every pair ties on score and on summed
  # count (1 + 3); rule (ii) chooses the pair whose child containing `a` is
  # smallest lexicographically: {a} < {a, b} < {a, c}.
  trees <- ape::read.tree(text = c("((a, b), c);", "((a, c), b);",
                                   "(a, (b, c));"))
  for (ord in list(1:3, c(3, 1, 2), c(2, 3, 1))) {
    expect_equal(cladeKeys(Hipstr(trees[ord])), c("a,b,c", "b,c"))
  }
})

test_that("Hipstr() and MrHipstr() attain the brute-force optimum", {
  for (cfg in list(c(4, 5, 3), c(5, 6, 4), c(5, 12, 6), c(6, 9, 5),
                   c(6, 30, 12))) {
    for (seed in 1:3) {
      trees <- randomSample(cfg[[1]], cfg[[2]], cfg[[3]], seed)
      tab <- hipTable(trees)
      hip <- Hipstr(trees)
      sc <- hipScore(hipTriples(hip), tab)
      expect_false(is.null(sc))
      expect_equal(sc[[2]], bruteOptimum(trees, FALSE)[[2]])
      expect_equal(attr(hip, "logCredibility"), sc[[2]])

      mr <- MrHipstr(trees)
      scMr <- hipScore(hipTriples(mr), tab)
      expect_equal(scMr, bruteOptimum(trees, TRUE))
      expect_equal(attr(mr, "logCredibility"), scMr[[2]])
    }
  }
})

test_that("HIPSTR scores at least as well as every input tree", {
  trees <- randomSample(14, 60, 15, 4321)
  tab <- hipTable(trees)
  inputs <- do.call(rbind, lapply(trees, function(tr)
    hipScore(hipTriples(tr), tab)))
  hip <- Hipstr(trees)
  mr <- MrHipstr(trees)
  expect_gte(attr(hip, "logCredibility"), max(inputs[, 2]) - 1e-12)
  scMr <- hipScore(hipTriples(mr), tab)
  best <- max(inputs[, 1])
  expect_gte(scMr[[1]], best)
  if (scMr[[1]] == best) {
    expect_gte(scMr[[2]], max(inputs[inputs[, 1] == best, 2]) - 1e-12)
  }
  # Clade credibilities agree with direct counts
  cred <- as.numeric(tab$count[hipTriples(hip)[, 1]]) / length(trees)
  expect_equal(sort(attr(hip, "credibility")), sort(cred))
})

test_that("result is independent of tree order and tip numbering", {
  trees <- randomSample(10, 40, 30, 99)
  ref <- Hipstr(trees)
  refMr <- MrHipstr(trees)
  set.seed(1)
  shuffled <- lapply(trees[sample(length(trees))], function(tr)
    TreeTools::RenumberTips(tr, sample(tr[["tip.label"]])))
  shuffled <- structure(shuffled, class = "multiPhylo")
  expect_equal(cladeKeys(Hipstr(shuffled)), cladeKeys(ref))
  expect_equal(cladeKeys(MrHipstr(shuffled)), cladeKeys(refMr))
  expect_equal(attr(Hipstr(shuffled), "logCredibility"),
               attr(ref, "logCredibility"))
})

test_that("Hipstr() is idempotent and handles degenerate input", {
  tree <- TreeTools::RandomTree(paste0("t", 1:9), root = TRUE)
  hip <- Hipstr(c(tree, tree, tree))
  expect_equal(cladeKeys(hip), cladeKeys(tree))
  expect_equal(attr(hip, "credibility"), rep(1, 8))
  expect_equal(attr(hip, "logCredibility"), 0)

  single <- Hipstr(tree)
  expect_equal(cladeKeys(single), cladeKeys(tree))
  expect_equal(cladeKeys(MrHipstr(list(tree))), cladeKeys(tree))

  twoLeaf <- ape::read.tree(text = "(a, b);")
  two <- Hipstr(list(twoLeaf, twoLeaf))
  expect_equal(two[["tip.label"]], c("a", "b"))
  expect_equal(two[["Nnode"]], 1L)
  expect_equal(attr(two, "credibility"), 1)

  oneLeaf <- ape::read.tree(text = "(a);")
  expect_identical(Hipstr(list(oneLeaf, oneLeaf)), oneLeaf)
  expect_null(Hipstr(list()))
  expect_equal(cladeKeys(Hipstr(list(tree, NULL, tree))), cladeKeys(tree))
})

test_that("output is a valid cladewise phylo with labels in input order", {
  trees <- randomSample(12, 25, 8, 7)
  hip <- Hipstr(trees)
  expect_equal(hip[["tip.label"]], TreeTools::TipLabels(trees[[1]]))
  expect_equal(nrow(hip[["edge"]]), 22L)
  expect_equal(hip[["edge"]][1, 1], 13L)
  expect_null(hip[["edge.length"]])
  expect_equal(ape::reorder.phylo(hip, "cladewise")[["edge"]], hip[["edge"]])
  expect_equal(TreeTools::NTip(hip), 12L)
  # Compressed tip labels
  compressed <- ape::.compressTipLabel(trees)
  expect_equal(cladeKeys(Hipstr(compressed)), cladeKeys(hip))
})

test_that("aliases are identical", {
  trees <- randomSample(7, 10, 4, 3)
  expect_identical(Hipster(trees), Hipstr(trees))
  expect_identical(MrHipster(trees), MrHipstr(trees))
})

test_that("Hipstr() rejects invalid input", {
  rooted <- ape::read.tree(text = "((a, b), (c, d));")
  expect_error(Hipstr(5), "list of trees")
  expect_error(Hipstr(list(rooted, "tree")), "list of trees")
  expect_error(Hipstr(list(rooted, ape::read.tree(text = "((a, b, c), d);"))),
               "fully resolved")
  expect_error(MrHipstr(list(rooted, ape::read.tree(text = "(a, b, (c, d));"))),
               "fully resolved")
  expect_error(Hipstr(list(rooted, ape::read.tree(text = "((a, b), (c, e));"))),
               "same leaves")
  expect_error(Hipstr(list(rooted, ape::read.tree(text = "((a, b), c);"))),
               "same leaves")
  expect_error(Hipstr(list(ape::read.tree(text = "((a, a), (c, d));"))),
               "unique")
  # Singleton (degree-two) node
  knuckle <- ape::read.tree(text = "(((a, b)), (c, d));")
  expect_error(Hipstr(list(rooted, knuckle)), "fully resolved")
})

test_that("C++ entry point validates its input", {
  e <- ape::read.tree(text = "((a, b), (c, d));")[["edge"]]
  idx <- 1:4
  expect_error(hipstrConsensus(list(), list(), 4L, FALSE), "No trees")
  expect_error(hipstrConsensus(list(e), list(), 4L, FALSE), "must match")
  expect_error(hipstrConsensus(list(e), list(idx), 1L, FALSE), "two tips")
  bad <- e
  bad[1, 1] <- 1L
  expect_error(hipstrConsensus(list(bad), list(idx), 4L, FALSE), "invalid edge")
  twoParents <- e
  twoParents[2, 2] <- twoParents[3, 2]
  expect_error(hipstrConsensus(list(twoParents), list(idx), 4L, FALSE),
               "more than one parent|fully resolved")
  expect_error(hipstrConsensus(list(e), list(c(1L, 2L, 3L, 9L)), 4L, FALSE),
               "invalid tip")
  expect_error(hipstrConsensus(list(rbind(c(5L, 1L), c(5L, 2L), c(6L, 3L),
                                          c(6L, 4L), c(7L, 5L), c(7L, 6L))),
                               list(idx), 4L, FALSE), NA)
  expect_error(hipstrConsensus(list(rbind(c(5L, 1L), c(5L, 6L), c(6L, 3L),
                                          c(6L, 7L), c(7L, 5L), c(7L, 2L))),
                               list(idx), 4L, FALSE), "no root")
  # Root plus a detached cycle 6 <-> 7
  expect_error(hipstrConsensus(list(rbind(c(5L, 1L), c(5L, 2L), c(6L, 7L),
                                          c(6L, 3L), c(7L, 6L), c(7L, 4L))),
                               list(idx), 4L, FALSE), "not connected")
})
