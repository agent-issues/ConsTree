# Compare Hipstr() / MrHipstr() with TreeAnnotator X (BEAST X 10.5.0) and
# sumt 4.1.5 on each data set written by make-datasets.R.
#
#   Rscript check-reference.R <lib> <dataDir> <outDir> <treeannotator> <sumt>
#
# For every data set and method, records: wall time of each tool; the log
# clade credibility (logCC) and number of majority clades of each tool's tree,
# all rescored by the independent referee (TreeTools::SplitFrequency on
# root-augmented trees); TreeAnnotator X's own report of its tree and its
# -target score of ConsTree's tree; rooted RF distances; and ConsTree's count
# of tie-break decisions along the returned tree.

args <- commandArgs(TRUE)
.libPaths(c(args[[1]], .libPaths()))
dataDir <- args[[2]]
outDir <- args[[3]]
taBin <- args[[4]]
sumtBin <- args[[5]]
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
source(file.path(dirname(sub("--file=", "", grep("--file=", commandArgs(FALSE),
                                                   value = TRUE))), "referee.R"))
library(ConsTree)
cat("ConsTree", format(packageVersion("ConsTree")), "\n")

Timed <- function(cmd, args, log) {
  t0 <- proc.time()[["elapsed"]]
  status <- system2(cmd, args, stdout = log, stderr = log)
  list(status = status, time = proc.time()[["elapsed"]] - t0,
       log = readLines(log, warn = FALSE))
}
Grab <- function(lines, pattern) {
  hit <- grep(pattern, lines, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(sub(".*:\\s*([-0-9.eE]+).*$", "\\1", hit[[length(hit)]]))
}

datasets <- read.csv(file.path(dataDir, "datasets.csv"))
rows <- list()
for (d in seq_len(nrow(datasets))) {
  name <- datasets$dataset[[d]]
  nex <- file.path(dataDir, paste0(name, ".nex"))
  trees <- structure(ReadTreesFile(nex), class = "multiPhylo")
  labels <- TreeTools::TipLabels(trees[[1]])
  forestA <- lapply(trees, .AddRootTip)
  for (method in c("hipstr", "mrhipstr")) {
    majority <- method == "mrhipstr"
    fn <- if (majority) MrHipstr else Hipstr
    t0 <- proc.time()[["elapsed"]]
    mine <- fn(trees)
    tMine <- proc.time()[["elapsed"]] - t0
    diag <- ConsTree:::hipstrConsensus(
      lapply(trees, `[[`, "edge"),
      lapply(trees, function(tr) match(tr[["tip.label"]], labels)),
      length(labels), majority)
    mineFile <- file.path(outDir, sprintf("%s_%s_constree.tre", name, method))
    plain <- mine
    plain[["node.label"]] <- NULL
    ape::write.tree(plain, mineFile)

    taOut <- file.path(outDir, sprintf("%s_%s_ta.nex", name, method))
    ta <- Timed(taBin, c("-type", method, "-heights", "mean", "-burninTrees",
                         "0", "-threads", "1", nex, taOut),
                file.path(outDir, sprintf("%s_%s_ta.log", name, method)))
    taTree <- ReadTreesFile(taOut)[[1]]
    taReported <- Grab(ta$log, "tree's log clade credibility")
    taBuild <- Grab(ta$log, "^\\[[0-9.]+ secs\\]")
    tgt <- Timed(taBin, c("-target", mineFile, "-heights", "mean",
                          "-burninTrees", "0", "-threads", "1", nex,
                          file.path(outDir, "target_out.nex")),
                 file.path(outDir, sprintf("%s_%s_target.log", name, method)))
    taScoresMine <- Grab(tgt$log, "Target tree's log clade credibility")

    sumtBase <- file.path(outDir, sprintf("%s_%s_sumt", name, method))
    st <- Timed("env", c("PYTHONHASHSEED=0", sumtBin,
                         if (majority) "--mrhip" else "--hip", "--noblen",
                         "-b", "0", "--cpus", "1", "--informat", "nexus",
                         "--outformat", "newick", "--nometa", "-n",
                         "--basename", sumtBase, nex),
                file.path(outDir, sprintf("%s_%s_sumt.log", name, method)))
    sumtTree <- ReadTreesFile(paste0(sumtBase, if (majority) ".mrhip" else ".hip"))
    if (inherits(sumtTree, "phylo")) sumtTree <- list(sumtTree)
    sumtTree <- sumtTree[[1]]

    sMine <- RefereeScore(mine, forestA)
    sTa <- RefereeScore(taTree, forestA)
    sSumt <- RefereeScore(sumtTree, forestA)
    rows[[length(rows) + 1L]] <- data.frame(
      dataset = name, method = method, nTip = length(labels),
      nTree = length(trees), nClade = diag$nClade, nPair = diag$nPair,
      logCC_constree = sMine$logCC, logCC_ta = sTa$logCC,
      logCC_sumt = sSumt$logCC,
      attr_logCC = attr(mine, "logCredibility"),
      ta_reported = taReported, ta_scores_constree = taScoresMine,
      nMaj_constree = sMine$nMajority, nMaj_ta = sTa$nMajority,
      nMaj_sumt = sSumt$nMajority,
      unobserved = sMine$nUnobserved + sTa$nUnobserved + sSumt$nUnobserved,
      cred_matches_referee = isTRUE(all.equal(sort(attr(mine, "credibility")),
                                              sMine$credibility)),
      rf_constree_ta = RootedRF(mine, taTree),
      rf_constree_sumt = RootedRF(mine, sumtTree),
      rf_ta_sumt = RootedRF(taTree, sumtTree),
      nTiedScore = diag$nTiedScore, nTiedCount = diag$nTiedCount,
      time_constree = tMine, time_ta = ta$time, time_ta_build = taBuild,
      time_sumt = st$time, status_ta = ta$status, status_sumt = st$status)
    print(rows[[length(rows)]][, c("dataset", "method", "logCC_constree",
                                   "logCC_ta", "logCC_sumt", "rf_constree_ta",
                                   "rf_constree_sumt", "nTiedScore",
                                   "time_constree", "time_ta", "time_sumt")])
    write.csv(do.call(rbind, rows), file.path(outDir, "reference-comparison.csv"),
              row.names = FALSE)
  }
}
