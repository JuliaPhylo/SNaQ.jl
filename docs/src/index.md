# SNaQ

[SNaQ](https://github.com/JuliaPhylo/SNaQ.jl) is a [Julia](http://julialang.org)
ements the SNaQ method by
[Solís-Lemus & Cécile Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)
to estimate a phylogenetic network from quartet concordance factors.
See the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package, which SNaQ depends on, for background on phylogenetic networks
and for how to get help.

## References

See them in
[bibtex format](https://github.com/juliaphylo/SNaQ.jl/blob/master/CITATION.bib).

For the SNaQ network inference method itself:
- Claudia Solís-Lemus and Cécile Ané (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  PLoS Genetics 12(3):e1005896. [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

## Manual

```@contents
Pages = ["man/installation.md",
         "man/ticr_howtogetQuartetCFs.md",
         "man/snaq_est.md",
         "man/fixednetworkoptim.md",
         "man/expectedCFs.md",
         "man/bootstrap.md",
         "man/parallelcomputation.md",
         "man/multiplealleles.md",
         "man/error_reporting.md"]
Depth = 1
```

For help on individual functions, see the library:

```@contents
Pages = map(file -> joinpath("lib", file), readdir("lib"))
Depth = 1
```
