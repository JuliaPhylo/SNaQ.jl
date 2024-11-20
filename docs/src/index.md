# SNaQ

[SNaQ](https://github.com/JuliaPhylo/SNaQ.jl) is a [Julia](http://julialang.org)
package that implements the SNaQ method by
[Solís-Lemus & Cécile Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)
to estimate a phylogenetic network from quartet concordance factors.
See the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package, which SNaQ depends on, for background on phylogenetic networks
and for how to get help.

## References

See them in
[bibtex format](https://github.com/juliaphylo/SNaQ.jl/blob/master/CITATION.bib).

for the SNaQ network inference method itself:
- Claudia Solís-Lemus and Cécile Ané (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  PLoS Genetics 12(3):e1005896. [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

## Manual

```@contents
Pages = [
    "man/installation.md",
    "man/inputdata.md",
    "man/ticr_howtogetQuartetCFs.md",
    "man/snaq_plot.md",
    "man/dist_reroot.md",
    "man/fixednetworkoptim.md",
    "man/expectedCFs.md",
    "man/bootstrap.md",
    "man/multiplealleles.md"
]
Depth = 3
```

For help on individual functions, see the library:

```@contents
Pages = [
    "lib/public.md",
    "lib/internal.md",
]
Depth = 3
```
