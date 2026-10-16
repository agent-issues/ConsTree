# HIPSTR / MrHIPSTR validation (`Hipstr()`, `MrHipstr()`)

Reference implementations:

- TreeAnnotator X, BEAST X v10.5.0: `dr.app.tools.treeannotator.HIPSTRTreeBuilder`
  (byte-identical on master as of 2026-10-09). Run as
  `treeannotator -type hipstr|mrhipstr -heights mean -burninTrees 0 -threads 1`.
- sumt 4.1.5 / phylotreelib 2.3.1 (`SummaryTreeBuilder.hipstr_tree`), run as
  `PYTHONHASHSEED=0 sumt --hip|--mrhip --noblen -b 0 --cpus 1`.

## Reference semantics (from the source) and ConsTree's choices

| | TreeAnnotator X | sumt | ConsTree |
|---|---|---|---|
| Majority threshold | `credibility > 0.5` | `freq > 0.5` | `2 * count > nTree` (same) |
| MrHIPSTR priority | `+1e10` per majority clade, in log space | `+1e10` | lexicographic (nMajority, logCC): exact limit, no rounding |
| Tie on score | first pair in Java `HashSet` order (its `breakTies` block never fires) | larger summed child frequency, then Python `set` order (hash-seed dependent) | larger summed child count, then lexicographic tip-set rule; tolerance 1e-12 relative |
| Polytomy in input | unsupported (Java `assert`, disabled by default: extra children silently ignored) | error | error |

## Scripts

- `brute-force.R <lib> <out.csv>`: exhaustive enumeration of all rooted binary
  trees (n = 3..8; 135 135 at n = 8) with an independent bitmask clade table.
- `make-datasets.R <lib> <ebov.trees> <outDir>`: EBOV (Dudas et al. 2017,
  1 610 tips x 1 000 trees, github.com/ebov/space-time
  `Analyses/Phylogenetic/Makona_1610_cds_ig.1K.trees.zip`), random-leaf prunings
  of it, and simulated posteriors.
- `check-reference.R <lib> <dataDir> <outDir> <treeannotator> <sumt>`: runs all
  three tools and rescores every output tree with the independent referee in
  `referee.R` (TreeTools::SplitFrequency on root-augmented trees).
- `results/`: CSVs from the 2026-10-09 run on Hamilton.

## Results (2026-10-09, ConsTree 1.0.1.9000)

- Brute force: 2 640 / 2 640 (1 320 data sets x 2 methods) returned an optimal
  tree inside the HIPSTR space; 687 had more than one optimal tree.
- Reference agreement: 14 data sets x 2 methods = 28 comparisons per reference.
  Log clade credibility identical (referee-rescored, difference 0) in 28/28 for
  both TreeAnnotator X and sumt; clade count identical to TreeAnnotator X's
  "Total unique clades" in 28/28; topology identical (rooted RF 0) in 23/28 vs
  TreeAnnotator X and 26/28 vs sumt. Every topological difference is between
  co-optimal trees (equal objective) at a tied choice, explained by the
  tie rules above.
- EBOV 1 610 x 1 000: HIPSTR logCC -3039.5956 (paper: -3039.60), 633/634
  majority clades; MrHIPSTR -3040.1036, 634/634.
