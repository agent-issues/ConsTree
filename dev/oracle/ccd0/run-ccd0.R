# Run ConsTree::CCD0() on a BEAST tree file; save tree, attributes and timing.
# usage: Rscript run-ccd0.R TREES OUT.rds
args <- commandArgs(TRUE)
source("dev/oracle/ccd0/beast-io.R"); source("dev/oracle/ccd0/referee.R")
library(ConsTree)
tRead <- system.time(trees <- ReadBeastTrees(args[1]))[["elapsed"]]
if (inherits(trees, "phylo")) trees <- c(trees)
labels <- TranslateOrder(args[1])
if (length(labels) != length(trees[[1]][["tip.label"]])) labels <- sort(trees[[1]][["tip.label"]])
# Index taxa in translate order, as BEAST does, so tie-breaks are comparable
trees <- ape::.uncompressTipLabel(trees)
trees[[1]] <- TreeTools::RenumberTips(trees[[1]], labels)
gc(reset = TRUE)
tRun <- system.time(map <- CCD0(trees))[["elapsed"]]
memR <- sum(gc()[, 6])   # max MB used by R since reset (Ncells+Vcells)
st <- readLines("/proc/self/status"); hwm <- as.numeric(gsub("\\D", "", st[grepl("^VmHWM", st)])) / 1024
saveRDS(list(tree = map, tRead = tRead, tRun = tRun, memR = memR, peakRSS = hwm,
             nTree = length(trees), nTip = length(labels)), args[2])
cat(sprintf("CCD0: n=%d k=%d clades=%d splits=%.0f read=%.2fs CCD0=%.2fs peakRSS=%.0fMB\n",
            length(labels), length(trees), attr(map, "nClades"), attr(map, "nSplits"), tRead, tRun, hwm))
