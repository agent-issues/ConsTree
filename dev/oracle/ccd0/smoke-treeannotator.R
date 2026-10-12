# Smoke test: TreeAnnotator 2.7.7 -topology CCD0 on the worked examples, and
# which taxon numbering its lexicographic tie-break uses.
args <- commandArgs(TRUE); ta <- args[1]; outDir <- args[2]
source("dev/oracle/ccd0/beast-io.R"); source("dev/oracle/ccd0/referee.R")
library(ConsTree)
run <- function(trees, labels, name) {
  f <- file.path(outDir, paste0(name, ".trees")); o <- file.path(outDir, paste0(name, ".tre"))
  WriteBeastNexus(trees, f, labels)
  st <- system2(ta, c("-topology", "CCD0", "-burnin", "0", f, o), stdout = TRUE, stderr = TRUE)
  ref <- ReadTreeAnnotator(o)
  ours <- CCD0(lapply(trees, TreeTools::RenumberTips, labels))
  cat(name, ": TreeAnnotator ", ape::write.tree(ref), " | ConsTree ", ape::write.tree(ours),
      " | RF ", RootedRF(ref, ours), "\n", sep = "")
}
fig1 <- ape::read.tree(text = c(rep("(((A:1,B:1):1,C:2):1,(D:1,E:1):2);", 3),
  rep("(((A:2,(B:1,C:1):1):1,D:3):1,E:4);", 2), rep("(((A:2,(B:1,C:1):1):1,E:3):1,D:4);", 2)))
run(fig1, LETTERS[1:5], "fig1")
fig2 <- ape::read.tree(text = c("(((A:1,B:1):1,C:2):1,D:3);", "(((C:1,D:1):1,A:2):1,B:3);"))
run(fig2, LETTERS[1:4], "fig2")
lex <- ape::read.tree(text = c("(A:3,((B:1,C:1):1,D:2):1);", "(((B:1,C:1):1,A:2):1,D:3);"))
run(lex, c("A", "B", "C", "D"), "lex_alpha")
run(lex, c("D", "B", "C", "A"), "lex_reversed")
run(lex, c("B", "D", "A", "C"), "lex_mixed")
relab <- ape::read.tree(text = c("(a:3,((c:1,d:1):1,B:2):1);", "(((c:1,d:1):1,a:2):1,B:3);"))
run(relab, c("a", "c", "d", "B"), "lex_case")
