# Frequency-difference consensus tree

`Frequency()` computes the frequency-difference consensus, which retains
each split that occurs more often than every split that conflicts with
it (Goloboff et al. 2003) .

## Usage

``` r
Frequency(trees)
```

## Arguments

- trees:

  A list of trees, or a `multiPhylo` object. All entries must share the
  same leaf labels.

## Value

`Frequency()` returns an object of class `phylo` denoting the frequency
consensus, rooted as in the first entry of `trees`.

## Details

The frequency-difference consensus is at least as resolved as the
majority-rule consensus
([`Majority()`](https://constree.github.io/reference/Majority.md)), and
is contained within the greedy consensus
([`Greedy()`](https://constree.github.io/reference/Greedy.md)).

This implementation builds on the FDCT algorithm of Jansson et al.
(2026) ; please cite that paper when using this method.

## References

Goloboff PA, Farris JS, Källersjö M, Oxelman B, Ramírez MJ, Szumik CA
(2003). “Improvements to resampling measures of group support.”
*Cladistics*, **19**(4), 324–332.
[doi:10.1111/j.1096-0031.2003.tb00376.x](https://doi.org/10.1111/j.1096-0031.2003.tb00376.x)
.  
  
Jansson J, Sung W, Tabatabaee SA, Yang Y (2026). “A faster algorithm for
constructing the frequency difference consensus tree.” *Journal of
Computer and System Sciences*, **161**, 103831.
[doi:10.1016/j.jcss.2026.103831](https://doi.org/10.1016/j.jcss.2026.103831)
.

## See also

Closely related:
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md).

Other consensus methods:
[`Adams()`](https://constree.github.io/reference/Adams.md),
[`Average()`](https://constree.github.io/reference/Average.md),
[`Greedy()`](https://constree.github.io/reference/Greedy.md),
[`Local()`](https://constree.github.io/reference/Local.md),
[`Loose()`](https://constree.github.io/reference/Loose.md),
[`Majority()`](https://constree.github.io/reference/Majority.md),
[`MajorityPlus()`](https://constree.github.io/reference/MajorityPlus.md),
[`Quartet()`](https://constree.github.io/reference/Quartet.md),
[`RStar()`](https://constree.github.io/reference/RStar.md),
[`Strict()`](https://constree.github.io/reference/Strict.md),
[`Transfer()`](https://constree.github.io/reference/Transfer.md)

## Examples

``` r
trees <- ape::as.phylo(0:5, 8)
Frequency(trees)
#> 
#> Phylogenetic tree with 8 tips and 6 internal nodes.
#> 
#> Tip labels:
#>   t2, t8, t7, t6, t3, t5, ...
#> 
#> Rooted; no branch length.
```
