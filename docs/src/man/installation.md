# Installation

## Install Julia

To install Julia, follow instructions [here](http://julialang.org/downloads/).

Julia code is compiled just-in-time, meaning that the first time you run a
function, it will be compiled at that moment. So, please be patient!
Future calls to the function will be much faster.
Trying out toy examples for the first calls is a good idea.

## Install SNaQ

To install the package, type inside Julia:
```julia
using Pkg
Pkg.add("SNaQ")
```
If you already installed the package and want the latest registered version,
do the following (which will update all of your packages):
```julia
Pkg.update()
```
It is important to update the package regularly as it is
undergoing constant development.

`Pkg.update()` will install the latest registered version, but there
could be other improvements in the `main` branch of the
repository. If you want to update to the latest unregistered version
of the package, you can do
`Pkg.add(PackageSpec(name="SNaQ", rev="main"))`
just beware that the latest changes could be not as robust.
If you want to go back to the registered package, you can do
`Pkg.free("SNaQ")`.

Similarly, you can pin a version of the package
`Pkg.pin("SNaQ")` so that `Pkg.update()` will not modify
it. You can always free a pinned package with
`Pkg.free("SNaQ")`. More on package management
[here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

The SNaQ package has dependencies like
[PhyloNetworks](https://github.com/JuliaPhylo/PhyloNetworks.jl)
(see the `Project.toml` file for the full list), but everything is installed automatically.

[PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl)
has utilities to visualize networks, and for interoperability,
such as to export networks to R (which can then be plotted via R).

To install:

```julia
using Pkg
Pkg.add("PhyloNetworks")
Pkg.add("PhyloPlots")
```

PhyloPlots also depends on PhyloNetworks, and has further dependencies
like
[RCall](https://github.com/JuliaInterop/RCall.jl)

To check that your installation worked, type this in Julia to load the package.
This is something to type every time you start a Julia session:

```@example install
using PhyloNetworks
using SNaQ
```

You can see a list of all the functions with
```julia
varinfo(SNaQ)
```
and press `?` inside Julia to switch to help mode,
followed by the name of a function (or type) to get more details about it.

## Test example

We show here small examples on how to get more
info on an object, what's its type, and how to manipulate objects.

For example, let's take an object `raxmlCF` created from reading in some data in the form of gene trees
(see more on the data in [Inputs for SNaQ](@ref)):

```@repl install
raxmltrees = joinpath(dirname(pathof(SNaQ)), "..","examples","raxmltrees.tre");
raxmlCF = readtrees2CF(raxmltrees);
```

Typing `varinfo()` will provide a list of objects and packages in memory,
including `raxmlCF` that we just created.
If we want to know the type of a particular object, we do:
```@repl install
typeof(raxmlCF)
```
which shows us that `raxmlCF` is of type [`DataCF`](@ref).
If we want to know about the attributes the object has, we can type `?` in Julia,
followed by `DataCF` for a description.
We can also ask for a list of all its attributes with

```@repl install
fieldnames(typeof(raxmlCF))
```
For example, we see that one attribute is `numQuartets`: its the number of 4-taxon subsets
in the data. To see what this number is:
```@repl install
raxmlCF.numQuartets
```
We also noticed an attribute `quartet`. It is a vector of Quartet objects inside `raxmlCF`, so
```@repl install
raxmlCF.quartet[2].taxon
```
will provide the list of taxon names for the second 4-taxon subset in the data.
To see the observed CF, we can type
```@repl install
raxmlCF.quartet[2].obsCF
```
We can verify the type with
```@repl install
typeof(raxmlCF.quartet[2])
```