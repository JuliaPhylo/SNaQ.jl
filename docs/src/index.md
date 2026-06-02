# SNaQ

[SNaQ](https://github.com/JuliaPhylo/SNaQ.jl) is a [Julia](http://julialang.org)
package as a part of the [JuliaPhylo](https://juliaphylo.github.io/JuliaPhyloWebsite/) software ecosystem.
 The package implements the SNaQ method by
[Solís-Lemus & Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)
to estimate a phylogenetic network from quartet concordance factors.
See the [PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
package, which SNaQ depends on, for background on phylogenetic networks
and for how to get help.

Join the PhyloNetworks google group for updates
[here](https://groups.google.com/g/juliaphylo-users).

More information about the pre-processing steps to get the input data for SNaQ can be found in [PhyloUtilities](https://juliaphylo.github.io/PhyloUtilities/) and the [snaq tutorial](https://solislemuslab.github.io/snaq-tutorial/).


## References

For the SNaQ network inference method itself:
- Claudia Solís-Lemus and Cécile Ané (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  PLoS Genetics 12(3):e1005896. [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

For improvements to SNaQ including multi-threading and weighted move sampling:
- Nathan Kolbow, Kevin Kong, Tyler Chafin, John Justison, Cécile Ané, and Claudia Solís-Lemus (2025).
  SNaQ.jl: Improved scalability for phylogenetic network inference.
  [doi:10.64898/2026.04.17.719296v2](https://doi.org/10.64898/2026.04.17.719296v2)

For inference with SNaQ beyond level-1 network space, gradient computations, and gradient-based optimization:
- Nathan Kolbow, Josh Justison, Claudia Solís-Lemus (2026).
  (Pending)

For PhyloNetworks package:
- Claudia Solís-Lemus, Paul Bastide and Cécile Ané (2017).
  PhyloNetworks: a package for phylogenetic networks. Molecular Biology and Evolution 34(12):3292–3298. [doi:10.1093/molbev/msx235](https://academic.oup.com/mbe/article/34/12/3292/4103410)


See the references in
[bibtex format](https://github.com/juliaphylo/SNaQ.jl/blob/master/CITATION.bib).

## Manual

```@contents
Pages =  [
"man/installation.md",
"man/snaq_est.md",
"man/network_restrictions.md",
"man/correlated_inheritance.md",
"man/fixednetworkoptim.md",
"man/expectedCFs.md",
"man/bootstrap.md",
"man/parallelcomputation.md",
"man/multiplealleles.md",
"man/error_reporting.md"
]
Depth = 3
```

## Library

For help on individual functions, see the library:

```@contents
Pages = map(file -> joinpath("lib", file), readdir("lib"))
Depth = 1
```
