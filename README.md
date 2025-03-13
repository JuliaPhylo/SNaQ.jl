# SNaQ: Maximum pseudolikelihood estimation of species network <img src="docs/src/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>

|**Documentation**| **Build Status** | **Code Coverage**   | **Style Guide** |
|:---------------:|:----------------:|:-------------------:|:----------------|
|[![stable][docs-stable-img]][docs-stable-url] [![dev][docs-dev-img]][docs-dev-url] | [![build][build-img]][build-url] [![PkgEval][pgkeval-img]][pgkeval-url] [![aqua][aqua-img]][aqua-url] | [![coverage][codecov-img]][codecov-url] | [![Code Style: Blue][style-img]][style-url] [![collaborative][colprac-img]][colprac-url]

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://JuliaPhylo.github.io/SNaQ.jl/stable/
[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://JuliaPhylo.github.io/SNaQ.jl/dev/

[build-img]: https://github.com/JuliaPhylo/SNaQ.jl/actions/workflows/CI.yml/badge.svg?branch=main
[build-url]: https://github.com/JuliaPhylo/SNaQ.jl/actions/workflows/CI.yml?query=branch%3Amain
[pgkeval-img]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/SNaQ.svg
[pgkeval-url]: https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/SNaQ.html
[aqua-img]: https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg
[aqua-url]: https://github.com/JuliaTesting/Aqua.jl

[codecov-img]: https://codecov.io/gh/JuliaPhylo/SNaQ.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaPhylo/SNaQ.jl

[style-img]: https://img.shields.io/badge/code%20style-blue-4495d1.svg
[style-url]: https://github.com/invenia/BlueStyle
<!-- ColPrac: Contributor's Guide on Collaborative Practices for Community Packages -->
[colprac-img]: https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet
[colprac-url]: https://github.com/SciML/ColPrac

## Overview

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and
inheritance probabilities, and a heuristic search in the space of
level-1 phylogenetic networks.

To get help, check

- the latest [SNaQ documentation](https://juliaphylo.github.io/SNaQ.jl/dev)
- the latest [PhyloNetworks documentation](https://juliaphylo.github.io/PhyloNetworks.jl/dev)
- [PhyloUtilities](https://juliaphylo.github.io/PhyloUtilities/) with a step-by-step tutorial from multi-locus sequences to necessary input for SNaQ
- [tutorial](https://solislemuslab.github.io/snaq-tutorial/) for snaq estimation (2023 workshop)
- the [JuliaPhylo google group](https://groups.google.com/forum/#!forum/juliaphylo-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.

## Citing

For the SNaQ method in particular, please cite
- Claudia Sol&iacute;s-Lemus and C&eacute;cile An&eacute; (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896.
  [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)

For PhyloNetworks package, please cite:
- Claudia Solís-Lemus, Paul Bastide and Cécile Ané (2017). 
  PhyloNetworks: a package for phylogenetic networks. Molecular Biology and Evolution 34(12):3292–3298. [doi:10.1093/molbev/msx235](https://academic.oup.com/mbe/article/34/12/3292/4103410)

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).

> [!NOTE]
> Much of this package was formerly part of PhyloNetworks v0.16.4 (and prior).
> PhyloNetworks v0.17, v1.0 (and later) will be stripped of functions implementing the SNaQ method.