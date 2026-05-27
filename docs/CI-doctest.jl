using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=pwd()))
Pkg.instantiate()

using Documenter: DocMeta, doctest
using SNaQ, PhyloNetworks
DocMeta.setdocmeta!(SNaQ, :DocTestSetup, :(using SNaQ, PhyloNetworks); recursive=true)
doctest(SNaQ)