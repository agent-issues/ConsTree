# Compare ConsTree::CCD0() with TreeAnnotator 2.7.7 (-topology CCD0) outputs.
# usage: Rscript compare.R RESULTDIR OUT.csv   (expects <id>.trees path list in
# RESULTDIR/manifest.tsv: id, trees file, ta tree file, ct rds, ta time file, ct time file)
args <- commandArgs(TRUE)
source("dev/oracle/ccd0/beast-io.R"); source("dev/oracle/ccd0/referee.R")
man <- read.delim(file.path(args[1], "manifest.tsv"), header = FALSE, stringsAsFactors = FALSE,
                  col.names = c("id", "trees", "ta", "ct", "taTime", "ctTime"))
TimeV <- function(f) {
  if (!file.exists(f)) return(c(NA, NA))
  x <- readLines(f, warn = FALSE)
  el <- sub(".*: ", "", x[grepl("Elapsed \\(wall clock\\)", x)])
  secs <- if (length(el)) sum(as.numeric(strsplit(el, ":")[[1]]) * 60^(rev(seq_along(strsplit(el, ":")[[1]])) - 1)) else NA
  rss <- as.numeric(sub(".*: ", "", x[grepl("Maximum resident set size", x)])) / 1024
  c(secs, if (length(rss)) rss else NA)
}
rows <- lapply(seq_len(nrow(man)), function(i) {
  m <- man[i, ]
  ct <- readRDS(m$ct)
  trees <- ReadBeastTrees(m$trees)
  k <- length(trees)
  freq <- CladeFreq(trees)
  ctTree <- Unquote(ct$tree)
  ta <- if (file.exists(m$ta) && file.info(m$ta)$size > 0) tryCatch(ReadTreeAnnotator(m$ta), error = function(e) NULL) else NULL
  lcCT <- LogCredibility(ctTree, freq, k)
  lcTA <- if (is.null(ta)) NA else LogCredibility(ta, freq, k)
  tt <- TimeV(m$taTime); tc <- TimeV(m$ctTime)
  data.frame(id = m$id, nTip = ct$nTip, nTree = k, nClades = attr(ctTree, "nClades"),
             nSplits = attr(ctTree, "nSplits"),
             rf = if (is.null(ta)) NA else RootedRF(ta, ctTree),
             logCredCT = lcCT, logCredTA = lcTA, diff = lcCT - lcTA,
             selfCheck = abs(lcCT - attr(ctTree, "logCredibility")),
             taBinary = if (is.null(ta)) NA else ape::is.binary(ta),
             ctCCD0s = ct$tRun, ctTotal_s = tc[1], ctRSS_MB = tc[2],
             taTotal_s = tt[1], taRSS_MB = tt[2])
})
res <- do.call(rbind, rows)
write.csv(res, args[2], row.names = FALSE)
print(res, digits = 6)
