# Scaling of ConsTree::CCD0() on diffuse synthetic samples (independent
# uniform random rooted trees: nearly every clade distinct -- a worst case for
# the number of clades m).  One configuration per call so peak RSS is clean.
# usage: Rscript scaling.R n k OUT.csv
args <- commandArgs(TRUE)
n <- as.integer(args[1]); k <- as.integer(args[2])
library(ConsTree)
set.seed(n * 7919 + k)
labs <- paste0("t", seq_len(n))
trees <- lapply(seq_len(k), function(i) ape::rtree(n, tip.label = labs, br = NULL))
st0 <- readLines("/proc/self/status"); base <- as.numeric(gsub("\\D", "", st0[grepl("^VmRSS", st0)])) / 1024
t1 <- system.time(map <- CCD0(trees, cladeProbability = FALSE))[["elapsed"]]
t2 <- system.time(map2 <- CCD0(trees))[["elapsed"]]
st <- readLines("/proc/self/status"); hwm <- as.numeric(gsub("\\D", "", st[grepl("^VmHWM", st)])) / 1024
row <- data.frame(n = n, k = k, nClades = attr(map, "nClades"), nSplits = attr(map, "nSplits"),
                  time_noProb_s = t1, time_withProb_s = t2, rssBefore_MB = base, peakRSS_MB = hwm)
print(row)
write.table(row, args[3], sep = ",", append = file.exists(args[3]),
            col.names = !file.exists(args[3]), row.names = FALSE)
