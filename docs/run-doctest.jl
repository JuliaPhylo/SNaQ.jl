# From the base directory, run with `julia --project=docs/ docs/run-doctest.jl`
using Pkg
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter: DocMeta, doctest
using SNaQ, PhyloNetworks
DocMeta.setdocmeta!(SNaQ, :DocTestSetup, :(using SNaQ, PhyloNetworks); recursive=true)
doctest(SNaQ)