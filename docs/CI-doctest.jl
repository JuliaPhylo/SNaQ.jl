using Pkg
Pkg.activate(@__DIR__)
Pkg.develop(PackageSpec(path=joinpath(@__DIR__, "..")))
Pkg.instantiate()

using Documenter: DocMeta, doctest
using SNaQ, PhyloNetworks
DocMeta.setdocmeta!(SNaQ, :DocTestSetup, :(using SNaQ, PhyloNetworks); recursive=true)
doctest(SNaQ)