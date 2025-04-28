# Deploy documentation locally

These notes need to be updated once all packages are registered.

Interactively in `docs`:

```julia
pkg> activate .
pkg> status # just to check
```
Note that `PhyloNetworks` is already pointing at the `dev` folder. This is because the snaq `Project.toml` lists `PhyloNetworks` from the repo dev branch.

```julia
pkg> dev https://github.com/JuliaPhylo/SNaQ.jl.git
```

We want `PhyloPlots` to plot things in the docs, but right now, `PhyloPlots` depends on `PhyloTraits`, so we need both.

In julia, but inside docs:

```julia
using Pkg
Pkg.add(PackageSpec(url="https://github.com/JuliaPhylo/PhyloTraits.jl", rev="dev01"))
Pkg.add(PackageSpec(name="PhyloPlots", rev="dev11"))
```

# Run doctests 

To run just the doctests, from the base directory:

```shell
julia --project=docs/ docs/run-doctest.jl
```

or if finer control is needed for dependencies, doctests can be run interactively while in the base directory:

```shell
pkg> activate ./docs
# do what you need to do to update dependencies / dev packages
julia> include("./docs/run-doctest.jl")
```