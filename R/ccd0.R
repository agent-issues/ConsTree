#' CCD0 maximum a posteriori summary tree
#'
#' `CCD0()` returns the maximum a posteriori (MAP) tree of the CCD0 tree
#' distribution \insertCite{Berling2025}{ConsTree}, a fully resolved point
#' estimate of a posterior sample of rooted trees.
#' `CCD0MAP()` is an alias.
#'
#' A conditional clade distribution (CCD) is a distribution over rooted
#' binary trees, each tree's probability being the product of the probabilities
#' of its clade splits \insertCite{Larget2013}{ConsTree}.
#' The CCD0 parameterization uses only the observed clade frequencies.
#' Let \eqn{f(C)} count the \eqn{k} input trees that contain clade \eqn{C}, and
#' let \eqn{\mathrm{cc}(C) = f(C) / k}{cc(C) = f(C) / k} be its clade
#' credibility.
#' The CCD0 contains every tree that can be assembled from observed clades:
#' a clade \eqn{C} may be split into any pair of observed clades
#' \eqn{A}, \eqn{B} with \eqn{A \cup B = C}{A + B = C}, whether or not
#' that split was itself observed.
#' (Restricting splits to those observed in the sample instead gives HIPSTR.)
#' Each such tree \eqn{T} has probability
#' \deqn{\Pr(T) = \frac{1}{Z}\prod_{C \in T} \mathrm{cc}(C)}{Pr(T) = prod cc(C) / Z}
#' where \eqn{Z} sums the product over all trees in the CCD0.
#'
#' The CCD0-MAP tree maximizes \eqn{\prod_{C \in T} \mathrm{cc}(C)}{prod cc(C)},
#' the criterion of the maximum clade credibility (MCC) tree, but over all
#' trees in the CCD0 rather than only the sampled trees.
#' It is found exactly by dynamic programming over clades in order of
#' increasing size.
#' Ties are broken as in the reference implementation (the 'CCD' package
#' for BEAST 2, used by TreeAnnotator `-topology CCD0`): the more
#' balanced split is preferred, then the split whose smaller child contains
#' the alphabetically first leaf that the two candidate children do not share.
#' As in BEAST, labels are compared character by character
#' (`sort(method = "radix")`), so the result does not depend on how the
#' input trees number their leaves.
#' Values within a relative tolerance of \eqn{10^{-12}} count as ties.
#'
#' Like any fully resolved summary tree, the CCD0-MAP tree can contain clades
#' with low support; a point estimate is not a statement of confidence.
#' Inspect the clade credibilities returned with the tree, or report a
#' majority-rule consensus ([`Majority()`]) alongside it.
#'
#' @section Computational cost:
#' Every split is visited by testing, for each observed clade \eqn{C}, the
#' observed clades that share the first leaf of \eqn{C}; the partner clade is
#' found by hashing.
#' Splits are not stored: memory is dominated by the clade table,
#' about \eqn{m (n / 8 + 40)}{m (n / 8 + 40)} bytes for \eqn{m} distinct clades
#' on \eqn{n} leaves, alongside the input trees themselves.
#' Running time grows with the number of candidate pairs, roughly
#' \eqn{O(m^2 / n)} for a diffuse sample and close to \eqn{O(m)} for a
#' concentrated one.
#'
#' @param trees A list of rooted trees of class `phylo`, or a `multiPhylo`
#' object, typically a posterior sample from which burn-in has been removed.
#' All trees must bear the same leaves.
#' Trees need not be binary: clades of polytomous trees are counted, but
#' the result is always fully resolved.
#' @param cladeProbability Logical; if `TRUE` (the default), compute the
#' probability of each clade of the MAP tree under the CCD0.
#' This roughly doubles the running time.
#'
#' @return `CCD0()` returns a rooted, fully resolved tree of class `phylo`,
#' without edge lengths.
#' Its `node.label` gives the clade credibility (sample frequency) of each
#' node's clade, as text.
#' Further attributes:
#' - `credibility`: numeric clade credibility of each internal node, in node
#'   order;
#' - `probability`: probability of each internal node's clade under the CCD0
#'   (`NA` if `cladeProbability = FALSE`);
#' - `logProbability`: natural logarithm of the probability of the tree under
#'   the CCD0;
#' - `logCredibility`: natural logarithm of the tree's clade credibility,
#'   \eqn{\sum \log \mathrm{cc}(C)}{sum log cc(C)};
#' - `logNTrees`: natural logarithm of the number of trees in the CCD0;
#' - `nClades`: number of distinct clades in the sample, including the leaves
#'   and the root;
#' - `nSplits`: number of clade splits in the CCD0.
#'
#' @examples
#' # A posterior sample whose most frequent tree is not the CCD0-MAP tree
#' # (Berling et al. 2025, fig. 1)
#' trees <- ape::read.tree(text = c(
#'   rep("(((A, B), C), (D, E));", 3),
#'   rep("(((A, (B, C)), D), E);", 2),
#'   rep("(((A, (B, C)), E), D);", 2)
#' ))
#' map <- CCD0(trees)
#' ape::write.tree(map)        # ((A,(B,C)),(D,E)): never sampled
#' attr(map, "credibility")
#' exp(attr(map, "logProbability"))  # 12 / 49
#'
#' @family consensus methods
#' @references \insertAllCited{}
#' @importFrom ape .uncompressTipLabel
#' @importFrom TreeTools RenumberTips TipLabels TreeIsRooted
#' @export
CCD0 <- function(trees, cladeProbability = TRUE) {
  if (inherits(trees, "phylo")) {
    trees <- list(trees)
  }
  if (!is.list(trees)) {
    stop("`trees` must be a list of trees or a `multiPhylo` object.")
  }
  if (!is.null(attr(trees, "TipLabel"))) {
    trees <- .uncompressTipLabel(trees)
  }
  trees <- c(trees)
  trees <- trees[!vapply(trees, is.null, logical(1))]
  nTree <- length(trees)
  if (nTree < 1L) {
    stop("`trees` contains no trees.")
  }
  if (!all(vapply(trees, inherits, logical(1), "phylo"))) {
    stop("Every element of `trees` must be a tree of class `phylo`.")
  }
  labels <- unname(TipLabels(trees[[1L]]))
  n <- length(labels)
  if (anyDuplicated(labels)) {
    stop("all tip labels must be unique")
  }
  if (any(vapply(trees[-1L], function(tr)
    !setequal(TipLabels(tr), labels) || length(TipLabels(tr)) != n,
    logical(1)))) {
    stop("all trees must have the same tip labels")
  }
  if (n < 2L) {
    stop("CCD0() needs trees with at least two leaves.")
  }
  rooted <- vapply(trees, TreeIsRooted, logical(1))
  if (!all(rooted)) {
    unrooted <- which(!rooted)
    stop("CCD0() summarizes rooted trees; tree ",
         paste(unrooted[seq_len(min(5L, length(unrooted)))], collapse = ", "),
         if (length(unrooted) > 5L) ", ...", " is unrooted.")
  }

  # Each tree on its OWN root, with tips 1..n in `labels` order.
  edgeList <- lapply(trees, function(tr) {
    if (!identical(tr[["tip.label"]], labels)) {
      tr[["tip.label"]] <- unname(tr[["tip.label"]])
      tr <- RenumberTips(tr, labels)
    }
    tr[["edge"]]
  })
  # Tie-breaking rank: ordinal (C-locale) label order, as BEAST's TreeParser
  tieRank <- match(labels, sort(labels, method = "radix")) - 1L
  res <- ccd0MapCpp(edgeList, n, tieRank, isTRUE(cladeProbability))

  credibility <- res[["count"]] / nTree
  tree <- structure(list(
    edge = res[["edge"]],
    tip.label = labels,
    Nnode = n - 1L,
    node.label = as.character(credibility)
  ), class = "phylo", order = "cladewise")
  attr(tree, "credibility") <- credibility
  attr(tree, "probability") <- res[["probability"]]
  attr(tree, "logProbability") <- res[["logCredibility"]] -
    res[["logNormalizer"]]
  attr(tree, "logCredibility") <- res[["logCredibility"]]
  attr(tree, "logNTrees") <- res[["logNTrees"]]
  attr(tree, "nClades") <- res[["nClades"]]
  attr(tree, "nSplits") <- res[["nSplits"]]
  # Return:
  tree
}

#' @rdname CCD0
#' @export
CCD0MAP <- CCD0
