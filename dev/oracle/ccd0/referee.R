# Independent referee for CCD0(): clade frequencies from ape::prop.part (ape's
# C code, not ConsTree's), keyed by sorted tip labels, so the CCD0 objective of
# ANY tree -- e.g. TreeAnnotator's -- can be scored without ConsTree.

.Keys <- function(pp, labels) {
  vapply(pp, function(x) paste(sort(labels[x]), collapse = "\r"), character(1))
}

CladeKeys <- function(tree) {
  .Keys(ape::prop.part(tree), tree[["tip.label"]])   # includes the root
}

CladeFreq <- function(trees) {
  trees <- ape::.compressTipLabel(trees)
  pp <- ape::prop.part(trees, check.labels = FALSE)
  structure(attr(pp, "number"), names = .Keys(pp, attr(pp, "labels")))
}

# Sum of log clade credibilities over the tree's clades (root contributes 0)
LogCredibility <- function(tree, freq, k) {
  f <- unname(freq[CladeKeys(tree)])
  if (anyNA(f)) return(-Inf)
  sum(log(f / k))
}

# Rooted RF distance: half the symmetric difference of clade sets
RootedRF <- function(t1, t2) {
  a <- CladeKeys(t1); b <- CladeKeys(t2)
  length(setdiff(union(a, b), intersect(a, b))) / 2
}

# Taxon order of a NEXUS translate block (BEAST numbering)
TranslateOrder <- function(file) {
  con <- file(file, "r"); on.exit(close(con))
  inT <- FALSE; out <- character(0)
  repeat {
    l <- readLines(con, n = 1L, warn = FALSE)
    if (!length(l)) break
    s <- trimws(l)
    if (grepl("^translate", s, ignore.case = TRUE)) { inT <- TRUE; next }
    if (inT) {
      if (grepl("^;", s)) break
      p <- strsplit(sub("[,;]$", "", s), "[[:space:]]+")[[1]]
      if (length(p) >= 2) out[as.integer(p[1])] <- gsub("^'|'$", "", paste(p[-1], collapse = " "))
      if (grepl(";$", s)) break
    }
    if (grepl("^tree ", s, ignore.case = TRUE)) break
  }
  out
}
