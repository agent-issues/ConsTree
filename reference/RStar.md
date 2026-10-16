# R\* consensus tree

`RStar()` returns the R\* consensus (Degnan et al. 2009) of a set of
rooted trees.

## Usage

``` r
RStar(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object. All entries must share the
  same leaf labels.

## Value

`RStar()` returns the consensus tree, an object of class `phylo`. It is
rooted by construction, but when the resolved triplets leave the deepest
divergence unresolved the root is a polytomy.

## Details

The R\* consensus is a rooted-triplet method. For every set of three
leaves it tallies, across the input trees, the three possible resolved
rooted triplets (`ab|c`, `ac|b`, `bc|a`) and keeps whichever appears
most frequently. Ties are not kept. The kept triplets form the set of
majority resolved triplets, \\R\_{maj}\\. Then \\R^\*\\ is the most
resolved tree that displays no resolved triplet outside \\R\_{maj}\\.

\\R^\*\\ is always a refinement of the majority-rule consensus: every
majority clade also appears in `RStar()`.

`RStar()` finds the strong clusters of \\R\_{maj}\\ among candidate
clusters built from triplet similarities (Jansson and Sung 2013; Jansson
et al. 2016) . For two trees it takes \\O(n^2)\\ time (Jansson et al.
2016) , testing each candidate with the criterion of Jansson and Sung
(2013) . For more trees it tallies every triplet directly (Bryant 2003)
, in \\O(kn^3)\\ time; this outpaces the \\O(n^2 \log^{k+2} n)\\
algorithm of Jansson et al. (2016) unless \\k\\ is very small. Set
`options(mc.cores = )` to count triplets on several threads.

## References

Bryant D (2003). “A classification of consensus methods for
phylogenetics.” In Janowitz MF, Lapointe F, McMorris FR, Mirkin B,
Roberts FS (eds.), *Bioconsensus*, volume 61 of *DIMACS Series in
Discrete Mathematics and Theoretical Computer Science*, 163–184.
American Mathematical Society.
[doi:10.1090/dimacs/061/11](https://doi.org/10.1090/dimacs/061/11) .  
  
Degnan JH, DeGiorgio M, Bryant D, Rosenberg NA (2009). “Properties of
consensus methods for inferring species trees from gene trees.”
*Systematic Biology*, **58**(1), 35–54.
[doi:10.1093/sysbio/syp008](https://doi.org/10.1093/sysbio/syp008) .  
  
Jansson J, Sung W (2013). “Constructing the R\* consensus tree of two
trees in subcubic time.” *Algorithmica*, **66**(2), 329–345.
[doi:10.1007/s00453-012-9639-1](https://doi.org/10.1007/s00453-012-9639-1)
.  
  
Jansson J, Sung W, Vu H, Yiu S (2016). “Faster algorithms for computing
the R\* consensus tree.” *Algorithmica*, **76**(4), 1224–1244.
[doi:10.1007/s00453-016-0122-2](https://doi.org/10.1007/s00453-016-0122-2)
.

## See also

Closely related:
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Local()`](https://constree.github.io/reference/Local.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Frequency()`](https://constree.github.io/reference/Frequency.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Transfer()`](https://constree.github.io/reference/Transfer.md)

## Examples

``` r
# Five trees whose majority signal recovers the species tree (((a,b),c),d):
trees <- c(
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, b), c), d);"),
  ape::read.tree(text = "(((a, c), b), d);"),
  ape::read.tree(text = "(((b, c), a), d);")
)

# (a, b) wins {a,b,c} by plurality (3 vs 1 vs 1)
ape::write.tree(RStar(trees))
#> [1] "(((a,b),c),d);"
```
