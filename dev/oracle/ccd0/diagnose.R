# Explain every topological disagreement between ConsTree::CCD0() and
# TreeAnnotator: for each clade present in both trees whose children differ,
# report both splits' subtree objectives, whether the products of integer clade
# counts are EXACTLY equal (prime factorization), and what the tie rule picks.
# usage: Rscript diagnose.R OUT.csv RESULTDIR [RESULTDIR ...]
args <- commandArgs(TRUE)
source("dev/oracle/ccd0/beast-io.R"); source("dev/oracle/ccd0/referee.R")
Factor <- function(x) {   # prime exponents of integer x, as named vector
  out <- integer(0); p <- 2L
  while (x > 1L) {
    while (x %% p == 0L) { nm <- as.character(p); out[nm] <- if (is.na(out[nm])) 1L else out[nm] + 1L; x <- x %/% p }
    p <- p + 1L
    if (p * p > x && x > 1L) { nm <- as.character(x); out[nm] <- if (is.na(out[nm])) 1L else out[nm] + 1L; break }
  }
  out
}
ProdExponents <- function(counts) {
  tot <- integer(0)
  for (f in counts) { e <- Factor(as.integer(f)); for (nm in names(e)) tot[nm] <- if (is.na(tot[nm])) e[[nm]] else tot[nm] + e[[nm]] }
  tot[order(as.integer(names(tot)))]
}
# Clades of tree as list: key -> children keys
Children <- function(tree) {
  tree <- TreeTools::Preorder(tree)
  n <- TreeTools::NTip(tree); tl <- tree[["tip.label"]]; edge <- tree[["edge"]]
  desc <- vector("list", n + tree[["Nnode"]]); for (i in seq_len(n)) desc[[i]] <- tl[i]
  for (e in rev(seq_len(nrow(edge)))) desc[[edge[e, 1]]] <- c(desc[[edge[e, 1]]], desc[[edge[e, 2]]])
  key <- vapply(desc, function(x) paste(sort(x), collapse = "\r"), "")
  kids <- split(edge[, 2], edge[, 1])
  res <- lapply(names(kids), function(p) sort(key[kids[[p]]]))
  names(res) <- key[as.integer(names(kids))]
  list(children = res, key = key)
}
SubtreeKeys <- function(ch, x) {  # internal clade keys in subtree rooted at x (incl. x)
  out <- x; st <- x
  while (length(st)) { y <- st[1]; st <- st[-1]; k <- ch[[y]]; if (!is.null(k)) for (z in k) if (!is.null(ch[[z]])) { out <- c(out, z); st <- c(st, z) } }
  out
}
rows <- list()
for (dir in args[-1]) {
  man <- read.delim(file.path(dir, "manifest.tsv"), header = FALSE, stringsAsFactors = FALSE,
                    col.names = c("id", "trees", "ta", "ct", "taTime", "ctTime"))
  for (i in seq_len(nrow(man))) {
    m <- man[i, ]
    if (!file.exists(m$ta) || file.info(m$ta)$size == 0) next
    ta <- ReadTreeAnnotator(m$ta); ct <- Unquote(readRDS(m$ct)$tree)
    if (RootedRF(ta, ct) == 0) next
    trees <- ReadBeastTrees(m$trees); k <- length(trees); freq <- CladeFreq(trees)
    a <- Children(ct)$children; b <- Children(ta)$children
    common <- intersect(names(a), names(b))
    diff <- common[vapply(common, function(x) !identical(a[[x]], b[[x]]), logical(1))]
    # top-most divergences only: not inside another divergent clade's subtree
    for (x in diff) {
      xs <- strsplit(x, "\r")[[1]]
      up <- vapply(setdiff(diff, x), function(y) all(xs %in% strsplit(y, "\r")[[1]]), logical(1))
      sa <- SubtreeKeys(a, x); sb <- SubtreeKeys(b, x)
      fa <- unname(freq[sa]); fb <- unname(freq[sb])
      objA <- sum(log(fa / k)); objB <- sum(log(fb / k))
      exact <- identical(ProdExponents(fa) * 1L, ProdExponents(fb) * 1L) && length(sa) == length(sb)
      sizes <- function(ch) sort(lengths(strsplit(ch, "\r")))
      small <- function(ch) { s <- lengths(strsplit(ch, "\r")); if (s[1] != s[2]) ch[which.min(s)] else sort(ch, method = "radix")[1] }
      rows[[length(rows) + 1L]] <- data.frame(id = m$id, cladeSize = lengths(strsplit(x, "\r")),
        nested = any(up), ctSplit = paste(sizes(a[[x]]), collapse = "|"), taSplit = paste(sizes(b[[x]]), collapse = "|"),
        objCT = objA, objTA = objB, objDiff = objA - objB, exactTie = exact,
        ctMoreBalanced = min(sizes(a[[x]])) > min(sizes(b[[x]])),
        sameBalance = min(sizes(a[[x]])) == min(sizes(b[[x]])))
    }
  }
}
res <- if (length(rows)) do.call(rbind, rows) else data.frame()
write.csv(res, args[1], row.names = FALSE)
print(res, digits = 15)
