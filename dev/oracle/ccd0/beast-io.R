# Helpers for validating ConsTree::CCD0() against TreeAnnotator 2 (-topology CCD0).
# Not shipped (dev/ is in .Rbuildignore).

# Write trees as a BEAST-2-style NEXUS file: taxa block, translate block with
# taxa numbered 1..n in `labels` order, one "tree STATE_i = ..." line per tree.
# Edge lengths are required by TreeAnnotator; missing ones are set to 1.
WriteBeastNexus <- function(trees, file, labels = TreeTools::TipLabels(trees[[1]])) {
  n <- length(labels)
  lines <- c("#NEXUS", "", "Begin taxa;", paste0("\tDimensions ntax=", n, ";"),
             "\t\tTaxlabels", paste0("\t\t\t", labels), "\t\t\t;", "End;",
             "Begin trees;", "\tTranslate",
             paste0("\t\t", seq_len(n), " ", labels,
                    ifelse(seq_len(n) < n, ",", "")), ";")
  out <- file(file, "w")
  on.exit(close(out))
  writeLines(lines, out)
  for (i in seq_along(trees)) {
    tr <- trees[[i]]
    if (is.null(tr[["edge.length"]])) tr[["edge.length"]] <- rep(1, nrow(tr[["edge"]]))
    tr[["tip.label"]] <- as.character(match(tr[["tip.label"]], labels))
    tr[["node.label"]] <- NULL
    writeLines(paste0("tree STATE_", i - 1L, " = ", ape::write.tree(tr)), out)
  }
  writeLines("End;", out)
}

# Temporary files go to $CCD0_TMP when set (compute-node /tmp can be small).
.Tmp <- function() tempfile(tmpdir = Sys.getenv("CCD0_TMP", tempdir()), fileext = ".nex")

# Strip NEXUS single quotes from tip labels, so that files whose translate
# block quotes taxon names compare equal to files that do not.
Unquote <- function(tr) {
  if (inherits(tr, "multiPhylo")) {
    tr <- ape::.uncompressTipLabel(tr)
    return(structure(lapply(tr, Unquote), class = "multiPhylo"))
  }
  tr[["tip.label"]] <- gsub("^'(.*)'$", "\\1", tr[["tip.label"]])
  tr
}

# Read a TreeAnnotator output tree, discarding [&...] annotations.
ReadTreeAnnotator <- function(file) {
  txt <- readLines(file, warn = FALSE)
  txt <- gsub("\\[[^]]*\\]", "", txt)
  tmp <- .Tmp()
  writeLines(txt, tmp)
  tr <- ape::read.nexus(tmp)
  unlink(tmp)
  if (inherits(tr, "multiPhylo")) tr <- tr[[1]]
  Unquote(tr)
}

# Read a BEAST posterior tree file (NEXUS, translate block) without
# annotations; returns multiPhylo with full labels.
ReadBeastTrees <- function(file) {
  txt <- readLines(file, warn = FALSE)
  txt <- gsub("\\[[^]]*\\]", "", txt)
  tmp <- .Tmp()
  writeLines(txt, tmp)
  on.exit(unlink(tmp))
  Unquote(ape::read.nexus(tmp))
}
