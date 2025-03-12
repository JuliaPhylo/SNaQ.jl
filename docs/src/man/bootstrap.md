```@setup bootstrap
using PhyloNetworks, SNaQ
mkpath("../assets/figures")
```
# Bootstrap

> This documentation pertains to SNaQ v1.0 as originally described in [Solís-Lemus & Cécile Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)

## Running a bootstrap analysis

There are two ways to do a bootstrap analysis.

First, we can use a table of quartet CFs estimated with credibility intervals,
such as if we used BUCKy.
The [TICR pipeline]() outputs a CF table with extra columns for credibility
intervals. We could then read that table and get bootstrap networks like this,
and tweak options as needed:

```julia
using DataFrames, CSV
df = CSV.read("tableCF_withCI.csv", DataFrame)
bootnet = bootsnaq(startnetwork, df, hmax=1, filename="bootstrap")
```

Alternatively, we can use bootstrap gene trees: one file of bootstrap trees per gene.
Here, the input is a text file that lists all the bootstrap files (one per gene).
We demonstrate this second option here.

The names of all our bootstrap files are listed in "BSlistfiles".
(ASTRAL can use the same file to do its own bootstrap, see
[PhyloUtilities](https://juliaphylo.github.io/PhyloUtilities/notebooks/Gene-Trees-RAxML.html)).

The function [`readmultinewick_files`](https://juliaphylo.github.io/PhyloNetworks.jl/stable/lib/public/#PhyloNetworks.readmultinewick_files-Tuple{AbstractString}) (from PhyloNetworks)
can read this list of file names, then
read each bootstrap file to get the bootstrap sample for each gene.
We can use them to sample input gene trees at random, one per gene,
and estimate a network from them. We ask the `bootsnaq` function
to repeat this resampling of bootstrap gene trees several times.

```julia
bootTrees = readmultinewick_files("BSlistfiles");
bootnet = bootsnaq(net0, bootTrees, hmax=1, nrep=10, runs=3,
                   filename="bootsnaq", seed=4321)
```

The bootstrap networks are saved in the `boostrap.out` file, so they
can be read in a new session with
`bootnet = readmultinewick("bootsnaq.out")`. To save the bootstrap networks to
a different file (perhaps after having re-rooted them with an
outgroup), we could do this: `writeMultiTopology(bootnet, "bootstrapNets.tre")`.

The example above asks for 10 bootstrap replicates,
which is definitely too few, to make the example run faster.
We might also increase the number of optimization runs (`runs`)
done for each bootstrap replicate. The following bootstrap was run with the
default 10 runs per replicate, and 100 bootstrap replicates,
and the 100 bootstrap networks come with the package:

```@example bootstrap
bootnet = readmultinewick(joinpath(dirname(pathof(SNaQ)), "..","examples","bootsnaq.out"));
length(bootnet)
```

If we used a specified list of quartets on the original data, we
should use that same list for the bootstrap analysis through the
option `quartetfile`.

## Support for tree and hybrid edges

Now that we have 100 bootstrap networks, we need to summarize
what they have in common (highly supported features) and what they
don't (areas of uncertainty).

We will use the PhyloNetworks functions to quantify [Network support](@extref).
