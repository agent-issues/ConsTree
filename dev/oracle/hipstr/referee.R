# Independent referee and I/O helpers for validating Hipstr() / MrHipstr()
# against TreeAnnotator X (BEAST X 10.5.0) and sumt 4.1.5.
# Not shipped (dev/ is .Rbuildignore'd).  Uses only ape and TreeTools.

# Read every tree in a NEXUS (optionally with a Translate block) or Newick
# file, discarding [&...] annotations.  Returns a list of phylo objects.
ReadTreesFile <- function(file) {
  txt <- readLines(file, warn = FALSE)
  txt <- gsub("\\[[^]]*\\]", "", txt)          # comments / annotations
  isNexus <- grepl("^#NEXUS", txt[[1]], ignore.case = TRUE)
  if (!isNexus) {
    return(ape::read.tree(text = paste(txt, collapse = "")))
  }
  whole <- paste(txt, collapse = "\n")
  trans <- NULL
  m <- regexpr("(?is)translate(.*?);", whole, perl = TRUE)
  if (m > 0) {
    block <- sub("(?is)^translate", "", regmatches(whole, m), perl = TRUE)
    block <- sub(";$", "", block)
    entries <- trimws(strsplit(block, "\n")[[1]])
    entries <- sub(",$", "", entries)
    entries <- entries[nzchar(entries) & entries != ";"]
    key <- sub("^(\\S+)\\s+.*$", "\\1", entries)
    val <- sub("^\\S+\\s+", "", entries)
    trans <- setNames(gsub("^'|'$", "", val), key)
  }
  treeLines <- grep("^\\s*tree\\s", txt, ignore.case = TRUE, value = TRUE)
  nwk <- sub("^[^=]*=\\s*", "", treeLines)
  trees <- lapply(nwk, function(x) ape::read.tree(text = x))
  if (!is.null(trans)) {
    trees <- lapply(trees, function(tr) {
      lab <- gsub("^'|'$", "", tr[["tip.label"]])
      tr[["tip.label"]] <- ifelse(lab %in% names(trans), trans[lab], lab)
      tr
    })
  }
  trees
}

# BEAST-style NEXUS with a Translate block, integer tip codes and [&R].
WriteBeastNexus <- function(trees, file) {
  labels <- trees[[1]][["tip.label"]]
  con <- file(file, "w")
  on.exit(close(con))
  writeLines(c("#NEXUS", "", "Begin taxa;",
               paste0("\tDimensions ntax=", length(labels), ";"),
               "\tTaxlabels", paste0("\t\t", labels), "\t\t;", "End;", "",
               "Begin trees;", "\tTranslate",
               paste0("\t\t", seq_along(labels), " ", labels,
                      c(rep(",", length(labels) - 1L), "")),
               "\t\t;"), con)
  for (i in seq_along(trees)) {
    tr <- trees[[i]]
    tr[["tip.label"]] <- as.character(match(tr[["tip.label"]], labels))
    tr[["node.label"]] <- NULL
    writeLines(paste0("tree STATE_", i - 1L, " = [&R] ",
                      ape::write.tree(tr)), con)
  }
  writeLines("End;", con)
}

.AddRootTip <- function(tr) {
  TreeTools::AddTip(tr, where = TreeTools::RootNode(tr), label = "ROOT__")
}

# Independent scorer: clade counts via TreeTools::SplitFrequency() on trees
# augmented with a tip attached at the root, so each unrooted split of the
# augmented tree is exactly one clade of the rooted tree.
#   forestA  -- lapply(trees, .AddRootTip), precomputed once per data set
# Returns logCC (sum of log credibilities over non-root clades; the root
# contributes log 1 = 0), nMajority (non-root clades with count > N / 2, plus 1
# for the root to match the C++ diagnostic), and the sorted credibilities.
RefereeScore <- function(tree, forestA) {
  N <- length(forestA)
  tr <- tree
  tr[["node.label"]] <- NULL
  attributes(tr) <- attributes(tr)[c("names", "class", "order")]
  cnt <- unname(TreeTools::SplitFrequency(.AddRootTip(tr), forestA))
  list(logCC = sum(log(cnt / N)),
       nMajority = sum(2 * cnt > N) + 1L,
       credibility = sort(c(1, cnt / N)),
       nUnobserved = sum(cnt == 0))
}

# Rooted Robinson-Foulds distance (count of clades in one tree but not the
# other, summed both ways).
RootedRF <- function(t1, t2) {
  a <- .AddRootTip(t1)
  b <- .AddRootTip(t2)
  shared <- sum(TreeTools::SplitFrequency(a, list(b)))
  nInternal <- TreeTools::NTip(t1) - 2L
  2L * (nInternal - shared)
}
