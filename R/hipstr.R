#' Highest independent posterior subtree reconstruction (HIPSTR)
#'
#' `Hipstr()` returns the HIPSTR summary tree
#' \insertCite{Baele2025}{ConsTree} of a sample of rooted trees -- typically a
#' Bayesian posterior sample of time trees, e.g. from BEAST.
#' `MrHipstr()` returns the majority-rule variant, MrHIPSTR.
#' `Hipster()` and `MrHipster()` are aliases.
#'
#' The _credibility_ of a clade is the proportion of input trees that contain
#' it, and the _clade credibility_ of a tree is the product of the
#' credibilities of its clades.
#' The maximum clade credibility (MCC) tree is the input tree with the highest
#' clade credibility.
#' HIPSTR searches a larger space: every fully resolved tree in which each
#' clade is split into a pair of child clades that were observed as siblings,
#' within that clade, in at least one input tree.
#' Each input tree belongs to this space, so the HIPSTR tree has a clade
#' credibility at least as high as the MCC tree, but it need not itself occur
#' in the sample.
#'
#' One pass over the input records every observed clade, its frequency, and
#' every observed pair of child clades.
#' A dynamic programme then visits clades in increasing size, giving each the
#' highest log clade credibility attainable by a subtree on that clade,
#' and the tree is read off the best pairs from the root down
#' \insertCite{Baele2025}{ConsTree}.
#' Clades are identified by a 128-bit hash of their leaves, so time and memory
#' are linear (expected) in the total size of the input.
#'
#' `MrHipstr()` maximises first the number of clades whose credibility
#' strictly exceeds one half, and only then the clade credibility.
#' This matches TreeAnnotator X and 'sumt', which add a reward of
#' \eqn{10^{10}}{1e10} to the log score of each such clade; `MrHipstr()` applies
#' the ordering exactly, avoiding the rounding error that the reward introduces
#' into sums of log credibilities.
#' The resulting tree usually contains every majority-rule clade, but this is
#' not guaranteed: a majority clade is included only if it can be reached
#' through observed pairs of child clades.
#' A clade found in exactly half the trees receives no priority.
#'
#' Pairs of child clades that tie on the objective (to within a relative
#' tolerance of \eqn{10^{-12}}{1e-12}, which absorbs rounding error in sums of
#' logarithms) are separated by preferring
#' (i) the pair whose two clades occur in more trees in total (the rule of
#' 'sumt'); then
#' (ii) the pair whose clade that contains the first tip (in the order of
#' `TipLabels(trees[[1]])`) has the lexicographically smallest sorted vector
#' of tip indices, a proper prefix counting as smaller.
#' The result is therefore deterministic and independent of the order of the
#' input trees.
#' TreeAnnotator X instead retains the first of the tied pairs in the
#' iteration order of a hash table, so it may return a different, but equally
#' optimal, tree.
#'
#' Input trees must be rooted and strictly bifurcating, as is every tree
#' sampled by BEAST; an error is raised otherwise.
#' A polytomy cannot be scored against the fully resolved alternatives
#' (it contains fewer clades, each of which reduces clade credibility),
#' and both reference implementations require bifurcating trees.
#' Unrooted trees should be rooted first, e.g. with [`TreeTools::RootTree()`].
#'
#' HIPSTR and MrHIPSTR trees are point estimates.
#' They are fully resolved even where the sample supports no resolution, so
#' many of their clades may have low credibility; unlike a majority-rule
#' consensus [`Majority()`], the topology alone does not convey uncertainty.
#' Inspect the credibility of each clade (see Value) before interpreting a
#' grouping.
#' Edge lengths are not estimated.
#'
#' @param trees A list of rooted, fully resolved trees of class `phylo`, or a
#' `multiPhylo` object, on a common set of leaves.
#'
#' @return `Hipstr()` and `MrHipstr()` return a rooted, fully resolved tree of
#' class `phylo`, without edge lengths, whose `node.label` gives the
#' credibility of the clade below each internal node.
#' The same values are stored, in order of node number (the root, with
#' credibility 1, first), in the numeric attribute `credibility`;
#' the attribute `logCredibility` holds the natural logarithm of the tree's
#' clade credibility (the sum of the logarithms of its clade credibilities).
#' A single input tree yields its own topology, with every credibility equal
#' to one; trees with fewer than two leaves are returned unchanged, and an
#' empty list returns `NULL`.
#'
#' @examples
#' trees <- ape::read.tree(text = c(
#'   "(((a, b), c), (d, e));",
#'   "(((a, b), c), (d, e));",
#'   "(((a, c), b), (d, e));",
#'   "(((a, b), (d, e)), c);",
#'   "((a, (b, c)), (d, e));"
#' ))
#' hip <- Hipstr(trees)
#' ape::write.tree(hip)
#' attr(hip, "credibility")
#' attr(hip, "logCredibility")
#'
#' # MrHIPSTR prioritizes clades found in more than half the trees
#' ape::write.tree(MrHipstr(trees))
#'
#' @family consensus methods
#' @seealso [`Majority()`] for the majority-rule consensus.
#' @references \insertAllCited{}
#' @importFrom ape .uncompressTipLabel
#' @importFrom TreeTools TipLabels
#' @export
Hipstr <- function(trees) {
  .Hipstr(trees, majority = FALSE)
}

#' @rdname Hipstr
#' @export
MrHipstr <- function(trees) {
  .Hipstr(trees, majority = TRUE)
}

#' @rdname Hipstr
#' @export
Hipster <- Hipstr

#' @rdname Hipstr
#' @export
MrHipster <- MrHipstr

.Hipstr <- function(trees, majority) {
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }
  if (!is.list(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  if (!is.null(attr(trees, "TipLabel"))) {
    trees <- .uncompressTipLabel(trees)
  }
  trees <- unclass(trees)
  trees <- trees[!vapply(trees, is.null, logical(1))]
  nTree <- length(trees)
  if (nTree == 0L) {
    return(NULL)
  }
  if (!all(vapply(trees, inherits, logical(1), "phylo"))) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  labels <- TipLabels(trees[[1]])
  n <- length(labels)
  if (anyDuplicated(labels)) {
    stop("All tip labels must be unique.")
  }
  tipIndex <- lapply(trees, function(tr) {
    idx <- match(tr[["tip.label"]], labels)
    if (length(idx) != n || anyNA(idx) || anyDuplicated(idx)) {
      stop("All `trees` must describe the same leaves.")
    }
    idx
  })
  if (n < 2L) {
    return(trees[[1]])
  }
  edges <- lapply(trees, `[[`, "edge")
  res <- hipstrConsensus(edges, tipIndex, n, majority)

  credibility <- res[["credibility"]]
  tree <- structure(list(
    edge = res[["edge"]],
    tip.label = labels,
    Nnode = n - 1L,
    node.label = as.character(credibility)
  ), class = "phylo", order = "cladewise")
  attr(tree, "credibility") <- credibility
  attr(tree, "logCredibility") <- res[["logCredibility"]]
  # Return:
  tree
}
