```@setup bootstrap
using PhyloNetworks, SNaQ
mkpath("../assets/figures")
```
# Bootstrap

## Running a bootstrap analysis

There are two ways to do a bootstrap analysis.

First, we can use a table of quartet CFs estimated with credibility intervals,
such as if we used BUCKy.
The [TICR pipeline](https://juliaphylo.github.io/PhyloUtilities/) outputs a CF table with extra columns for credibility
intervals. We could then read that table and get bootstrap networks like this,
and tweak options as needed:

```julia
using DataFrames, CSV
df = CSV.read("tableCF_withCI.csv", DataFrame)
bootnet = bootsnaq(startnetwork, df, hmax=1, filename="bootstrap")
```

However, it is not necessary to use BUCKy however and a posterior sample of gene trees to calculate confidence intervals on quartet concordance factors. The function `confintqCF_genetrees` allows to use both a collection of gene tree point estimates (e.g., maximum likelihood gene trees), and then a sample of bootstrap replicates, one for each gene tree. A simple way to estimate both is using IQTREE. Imagine we have 10 gene trees, and let's suppose that we have all of our maximum likelihood gene trees in a file named `gts.tre`. Also, we have ten files, one per gene tree (named `gt1.ufboot`, `gt2.ufboot`, etc., all in a directory `bootstrap_gts`), containing 1000 ultrafast bootstrap replicates. Let us assume we have the following directory structure:

```bash
.
├── maximum_likelihood
│   └── gts.tre
└── bootstrap_gts
    ├── gt01.ufboot
    ├── gt02.ufboot
    ├── gt03.ufboot
    ├── gt04.ufboot
    ├── gt05.ufboot
    ├── gt06.ufboot
    ├── gt07.ufboot
    ├── gt08.ufboot
    ├── gt09.ufboot
    └── gt10.ufboot

```

We can calculate the quartet concordance factor point estimates and their confidence intervals atthe 0.95 level:

```julia
using SNaQ, PhyloNetworks, PhyloPlots, DataFrames, CSV

# read the gene tree mles
gt_mles = readmultinewick("maximum_likelihood/gts.tre")

# read the gene tree bootstrap replicates, one file per gene tree
gt_bootstrap_path = filter(contains("ufboot"), readdir("bootstrap_gts", join=true))
gt_bootstrap = map(x -> readmultinewick(x), gt_bootstrap_path)

qcf_ci = confintqCF_genetrees(gt_mles, gt_bootstrap, 0.95, true)
15×14 DataFrame
 Row │ t1      t2      t3      t4      CF12_34  CF12_34_lo  CF12_34_hi  CF13_24  CF13_24_lo  CF13_24_hi  CF14_23  CF14_23_lo  CF14_23_hi  ngenes
     │ String  String  String  String  Float64  Float64     Float64     Float64  Float64     Float64     Float64  Float64     Float64     Float64
─────┼────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
   1 │ A       B       C       D           0.6      0.4         0.8         0.2         0.0         0.3      0.2         0.1         0.5     10.0
   2 │ A       B       C       E           0.7      0.4         0.8025      0.1         0.0         0.3      0.2         0.0         0.5     10.0
   3 │ A       B       D       E           0.6      0.5         0.8         0.2         0.1         0.3      0.2         0.0         0.3     10.0
   4 │ A       C       D       E           0.5      0.3         0.7         0.4         0.2         0.5      0.1         0.0         0.3     10.0
   5 │ B       C       D       E           0.6      0.4         0.8         0.3         0.1         0.4      0.1         0.0         0.3     10.0
   6 │ A       B       C       O           0.4      0.2         0.7         0.2         0.0         0.4      0.4         0.1         0.6     10.0
   7 │ A       B       D       O           0.4      0.2         0.7         0.2         0.1         0.4      0.4         0.1         0.6     10.0
   8 │ A       C       D       O           0.2      0.1         0.6         0.3         0.1         0.5      0.5         0.2         0.6     10.0
   9 │ B       C       D       O           0.4      0.2         0.7         0.2         0.0         0.4      0.4         0.1         0.6     10.0
  10 │ A       B       E       O           0.4      0.2         0.7         0.2         0.0         0.4      0.4         0.1         0.6     10.0
  11 │ A       C       E       O           0.2      0.0975      0.5         0.2         0.0         0.4      0.6         0.3         0.7     10.0
  12 │ B       C       E       O           0.3      0.1         0.6         0.2         0.0         0.4      0.5         0.2         0.7     10.0
  13 │ A       D       E       O           0.2      0.1         0.5         0.2         0.0         0.4      0.6         0.3         0.7     10.0
  14 │ B       D       E       O           0.1      0.0         0.4         0.4         0.1         0.5      0.5         0.3         0.8     10.0
  15 │ C       D       E       O           0.1      0.0         0.4         0.5         0.1         0.6      0.4         0.2         0.7     10.0
```

The `DataFrame` `qcf_ci` can then be used as the input for the `bootsnaq` function:

```julia
# read in the network to bootstrap
net = readnewick("network.out")
bootsnaq(net, qcf_ci, hmax=1, filename="bootstrap", runs=10, nrep=10)
```

Please note that it is not necessary to have all the point estimates in the same file, as `readmultinewick` also
accepts a vector with a list of files as input.
    
A key difference between using a vector of gene trees (see below) and the table with
confidence intervals is that such table does not assume complete taxon sampling across gene trees
and therefore can incorporate incomplete taxon sampling, that is, not all of the species
need to be in all gene trees. This is a key feature as this sort of missing data is quite
common in empirical datasets. If your dataset is like this, and not all species are in all gene trees,
it is necessary to use the table as described here. On the other hand, if there is complete species sampling
across gene trees, you may either use the table or a sample of trees as described below.
    
Alternatively, we can use bootstrap gene trees: one file of bootstrap trees per gene.
Here, the input is a text file that lists all the bootstrap files (one per gene).
Please note that when using bootstrap trees themselves, these should have all the species
in all the trees, if there are missing species the bootstrapped network
will only include those present in all the gene trees.
We demonstrate this second option here.

The names of all our bootstrap files are listed in "BSlistfiles".
(ASTRAL can use the same file to do its own bootstrap, see
[PhyloUtilities](https://juliaphylo.github.io/PhyloUtilities/notebooks/Gene-Trees-RAxML.html)).

The `PhyloNetworks` function [`PhyloNetworks.readmultinewick_files`](@extref)
can read this list of file names, then
read each bootstrap file to get the bootstrap sample for each gene.
We can use them to sample input gene trees at random, one per gene,
and estimate a network from them. We ask the [`bootsnaq`](@ref) function
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

The example above creates 10 bootstrap replicates,
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

We can use PhyloNetworks functions to quantify [Network support](@extref).
