# SNaQ: Maximum pseudolikelihood estimation of species network <img src="docs/src/snaq.png" align=right title="SNaQ logo" width=262.5 height=111>

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JuliaPhylo.github.io/SNaQ.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JuliaPhylo.github.io/SNaQ.jl/dev/)
[![Build Status](https://github.com/JuliaPhylo/SNaQ.jl/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/JuliaPhylo/SNaQ.jl/actions/workflows/CI.yml?query=branch%3Amaster)
[![Coverage](https://codecov.io/gh/JuliaPhylo/SNaQ.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaPhylo/SNaQ.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![PkgEval](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/SNaQ.svg)](https://JuliaCI.github.io/NanosoldierReports/pkgeval_badges/S/SNaQ.html)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

## Overview

SNaQ implements the statistical inference method in Sol&iacute;s-Lemus and An&eacute;
[(2016)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and
inheritance probabilities, and a heuristic search in the space of
level-1 phylogenetic networks.

To get help, check

- the latest [SNaQ documentation](https://juliaphylo.github.io/SNaQ.jl/dev)
- the latest [PhyloNetworks documentation](https://juliaphylo.github.io/PhyloNetworks.jl/dev)
- a [wiki](https://github.com/juliaphylo/PhyloNetworks.jl/wiki) for a
  step-by-step tutorial with background on networks (last revised 2022)
- [tutorial](https://cecileane.github.io/networkPCM-workshop/) for
  comparative methods, including network calibration (2023 workshop)
- the [google group](https://groups.google.com/forum/#!forum/phylonetworks-users)
  for common questions. Join the group to post/email your questions,
  or to receive information on new versions, bugs fixed, etc.

## Citing

See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
For the SNaQ method in particular, please cite

- Claudia Sol&iacute;s-Lemus and C&eacute;cile An&eacute; (2016).
  Inferring Phylogenetic Networks with Maximum Pseudolikelihood under Incomplete Lineage Sorting.
  [PLoS Genet](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896)
  12(3):e1005896.
  [doi:10.1371/journal.pgen.1005896](https://doi.org/10.1371/journal.pgen.1005896)
