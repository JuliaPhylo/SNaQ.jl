```@setup snaqplot
using PhyloNetworks, SNaQ
mkpath("../assets/figures")
exampledir = joinpath(dirname(pathof(SNaQ)), "..","examples")
raxmltrees = joinpath(exampledir,"raxmltrees.tre")
raxmlCF = readtrees2CF(raxmltrees, writeTab=false, writeSummary=false)
astralfile = joinpath(exampledir,"astral.tre")
astraltree = readmultinewick(astralfile)[102] # 102th tree = last tree here
net0 = readnewick(joinpath(exampledir,"net0.out"))
net1 = readnewick(joinpath(exampledir,"net1.out"))
rotate!(net1, -6)
net2 = readnewick(joinpath(exampledir,"net2.out"))
net3 = readnewick(joinpath(exampledir,"net3.out"))
loglik!(net0, 53.53150526187732)
loglik!(net1, 28.31506721890958)
loglik!(net2, 28.31506721890957)
loglik!(net3, 28.315067218909626)
```
# Network estimation

SNaQ implements the statistical inference method in
[Solís-Lemus & Ané 2016](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005896).
The procedure involves a numerical optimization of branch lengths and inheritance
probabilities and a heuristic search in the space of phylogenetic networks.

We suggest that you create a special directory for running these examples,
where input files can be downloaded and where output files will be
created (with estimated networks for instance). Enter this directory
and run Julia from there.

## Inputs for SNaQ 

SNaQ uses has two main inputs for estimating phylogenetic networks: concordance factors (CFs) and a starting tree (or network).

### Concordance factors

Concordance factors denote frequency of each quartet topology 
present among the gene trees, which can be estimated using
[MrBayes](http://mrbayes.sourceforge.net) or
[RAxML](http://sco.h-its.org/exelixis/software.html) for example. 

This [pipeline](https://juliaphylo.github.io/PhyloUtilities/) can be used to obtain the table of
quartet CF needed as input for SNaQ
(see also the [snaq tutorial](https://solislemuslab.github.io/snaq-tutorial/)).
It starts from the sequence alignments,
runs MrBayes and then BUCKy (both parallelized), producing the
table of estimated CFs and their credibility intervals.

#### CFs from gene trees 

Suppose you have a file with a list of gene trees in parenthetical
format called `raxmltrees.tre`.
You can access the example file of input trees
[here](https://github.com/juliaphylo/SNaQ/blob/main/examples/raxmltrees.tre)
(or
[here](https://raw.githubusercontent.com/juliaphylo/SNaQ/main/examples/raxmltrees.tre)
for easier download).

Do not copy-paste into a "smart" text-editor. Instead, save the file
directly into your working directory using "save link as" or "download linked file as".
This file contains 30 gene trees, each in parenthetical format on 6 taxa
like this (with rounded branch lengths):

`(E:0.038,((A:0.014,B:0.010):0.010,(C:0.008,D:0.002):0.010):0.025,O:0.078);`

If `raxmltrees.tre` is in your working directory, you can
read in all gene trees and directly summarize them by a list
of quartet CFs (proportion of input trees with a given quartet):
```@repl qcf
using PhyloNetworks,SNaQ
raxmltrees = joinpath(dirname(pathof(SNaQ)), "..","examples","raxmltrees.tre");
raxmlCF = readtrees2CF(raxmltrees) # read in the file and produce a "DataCF" object

```

In this table (`tableCF.txt`), each 4-taxon set is listed in one row.
The 3 "CF" columns gives the proportion of genes that has
each of the 3 possible trees on these 4 taxa.

When there are many more taxa, the number of quartets
might be very large and we might want to use a subset to speed things up.
Here, if we wanted to use a random sample of 10 quartets
instead of all quartets, we could do:

`raxmlCF = readtrees2CF(raxmltrees, whichQ="rand", numQ=10, CFfile="tableCF10.txt")`

Be careful to use a `numQ` value smaller than the total number of possible
4-taxon subsets, which is *n choose 4* on *n* taxa (e.g. 15 on 6 taxa).
To get a predictable random sample, you may set the seed with
`using Random; Random.seed!(12321)`
(for instance) prior to sampling the quartets as above.
Additionally, providing a file name for the optional argument `CFile` saves the quartet CFs to file for later use.

#### CFs from large datasets
When we want to get *all* quartet CFs, 
the `readtrees2CF` is *much* slower than the PhyloNetworks function [`countquartetsintrees`](https://juliaphylo.github.io/PhyloNetworks.jl/stable/lib/public/#PhyloNetworks.countquartetsintrees)
to read in trees and calculate the quartet CFs observed in the trees:

```@repl qcf
trees = readmultinewick(raxmltrees)
q,t = countquartetsintrees(trees);
nt = tablequartetCF(q,t); # named tuple
using DataFrames
df = DataFrame(nt, copycols=false); # convert to a data frame, without copying the column data
show(df, allcols=true) # data frames are displayed much more nicely than named tuples
```


#### Reading CFs directly from file 

If we already have a table of quartet concordance factors (CFs)
saved as a table in this format

| Taxon1 | Taxon2 | Taxon3 | Taxon4 | CF12_34 | CF13_24 | CF14_23
|:-------|:-------|:-------|:-------|:--------|:--------|:-------
| D      | A| E | O|   0.565 |       0.0903 |       0.3447
| ...    |  |   |  |         |              |       ...

then we could read it in one step using the`readtableCF` function.

Concordance factors (CF), i.e. gene tree frequencies, for each
4-taxon subset can be obtained from [BUCKy](http://www.stat.wisc.edu/~ane/bucky/)
to account for gene tree uncertainty.

An example file comes with the package, available
[here](https://github.com/juliaphylo/SNaQ/blob/main/examples/buckyCF.csv)
or
[here](https://raw.githubusercontent.com/juliaphylo/SNaQ/main/examples/buckyCF.csv).

```@repl qcf
buckyCFfile = joinpath(dirname(pathof(SNaQ)), "..","examples","buckyCF.csv");
buckyCF = readtableCF(buckyCFfile)
```
The same thing could be done in 2 steps:
first to read the file and convert it to a 'DataFrame' object,
and then to convert this DataFrame into a DataCF object.
```@repl qcf
using CSV, DataFrames
dat = CSV.read(buckyCFfile, DataFrame);
first(dat, 6) # to see the first 6 rows
buckyCF = readtableCF(dat)
tablequartetCF(buckyCF)
```
In the input file, columns need to be in the right order:
with the first 4 columns giving the names of the taxa in each 4-taxon set.
The CF values are assumed to be in columns named "CF12_34", etc.,
or else in columns 5,6,7.
If available, a column named "ngenes" will be taken to have the
the number of genes for each 4-taxon subset.

### Starting tree

The other input for SNaQ is a starting tree (or network) to be used as a starting point in optimization.


If we have a tree for the data set at hand,
it can be used as a starting point for the optimization.
From our gene trees, we estimated a species tree with
[ASTRAL](https://github.com/smirarab/ASTRAL/blob/master/astral-tutorial.md).
This tree comes with the package in file `astral.tre`
[here](https://github.com/juliaphylo/SNaQ/blob/main/examples/astral.tre).
This file has 102 trees: 100 bootstrap species trees,
followed by their greedy consensus,
followed by the best tree on the original data.
It's this last tree that we are most interested in.
We can read it with
```@example qcf
astralfile = joinpath(dirname(pathof(SNaQ)), "..","examples","astral.tre");
astraltree = readmultinewick(astralfile)[102] # 102th tree: last tree here
```

Instead of a starting tree (`astraltree` in this case), we can start the optimization 
in SNaQ with a network, but this network needs to be "level-1".
Note that all trees and all networks with 1 hybridization are of level 1.
To make sure that a network with 2 or more hybridizations is of level 1,
we can read it in with
`readnewicklevel1` (which also unroots the tree, resolves polytomies,
replaces missing branch lengths by 1 for starting values etc.):
```julia
T=readnewicklevel1("startNetwork.txt")
```
Here `startNetwork.txt` is a hypothetical file: replace this by
the name of a file that contains your network of interest.

## Estimating a network

After we have [Inputs for SNaQ](@ref), we can estimate the network using the
input data `raxmlCF` and starting from tree (or network) `astraltree`.
We first set `hmax=0` to impose the constraint of at most 0 hybrid node,
that is, we estimate a tree.
```julia
net0 = snaq!(astraltree,raxmlCF, hmax=0, filename="net0", seed=1234)
```
Part of the screen output shows this:

    MaxNet is (C,D,((B,A):1.395762055180493,(O,E):0.48453400554506426):10.0);
    with -loglik 53.53150526187732

This parenthetical (extended Newick) description is not very
human-friendly, so we use [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl) to plot the tree
(more about plotting networks can be found in the PhyloNetworks guide [Network Visualization](@extref) 
and in [PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl) documentation).

```@example snaqplot
using PhyloPlots
using RCall # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('snaqplot_net0_1.svg'), width=4, height=3)" # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net0);
R"dev.off()"; # hide
nothing # hide
```
![net0_1](../assets/figures/snaqplot_net0_1.svg)

We can use this tree as a starting point to search for the best
network allowing for at most `hmax=1` hybridizations (which is the default if we do
not specify a `hmax` value).
```julia
net1 = snaq!(net0, raxmlCF, hmax=1, filename="net1", seed=2345)
```
part of screen output:

    best network and networks with different hybrid/gene flow directions printed to .networks file
    MaxNet is (C,D,((O,(E,#H7:::0.19558838614943078):0.31352437658618976):0.6640664399202987,(B,(A)#H7:::0.8044116138505693):10.0):10.0);
    with -loglik 28.31506721890958

We can visualize the estimated network and its inheritance values γ, which
measure the proportion of genes inherited via each parent at a reticulation event
(e.g. proportion of genes inherited via gene flow).
```@example snaqplot
R"svg(name('snaqplot_net1_1.svg'), width=4, height=3)"; # hide
R"par"(mar=[0,0,0,0]); # hide
plot(net1, showgamma=true);
R"dev.off()"; # hide
nothing # hide
```
![net1_1](../assets/figures/snaqplot_net1_1.svg)

This network has A as a hybrid, 80.4% sister to B,
and 19.6% sister to E (which is otherwise sister to O).
C & D are sister to each other.

Note that SNaQ infers an unrooted semi-directed network; 
the lack of rooting is depicted visually with a trifurcationon on the leftmost side of the plot.
The direction of hybrid edges can be inferred,
but the direction of tree edges cannot be inferred.
To obtain a representative visualization,
it is best to root the network first, using one or more outgroup.
PhyloNetworks has a guide on [Re-rooting trees and networks](@extref) for this.
If your outgroup conflicts with the direction of reticulations
in the estimated network, see the section
[Candidate networks compatible with a known outgroup](@ref).

We can also check the output files created by `snaq!`:
```julia
less("net1.err") # would provide info about errors, if any
less("net1.out") # main output file with the estimated network from each run
less("net1.networks") # extra info
```


when viewing these result files with `less`
within Julia, use arrows to scroll down and type `q` to quit viewing the files.
- The `.networks` file contains a list of networks that are slight modifications
of the best (estimated) network `net1`. The modifications change the direction
of one reticulation at a time, by moving the placement of one hybrid node to another
node inside the same cycle.
For each modified network, the pseudolikelihood score was calculated
(the `loglik` or `-Ploglik` values give a pseudo deviance actually).

- The `.out` file contains the best network among all runs, and the best
  network per run, includes also the pseudolikelihood score and the
  computation time.

- The `.log` file contains a description of each run, convergence criterion, and seed information.

- The `.err` file has seed information on runs that failed, empty when nothing failed. In the case of a failed run, 
you could run the `snaqDebug` function on the same settings that caused the error (to help us debug).


The function name `snaq!` ends with ! because it modifies the argument `raxmlCF`
by including the expected CF. Type `?` then `snaq!` to get help on that function.

The main output file, here `net1.out` (or `snaq.out` by default) has the estimated
network in parenthetical format, but we can also print it directly to the screen:
```@repl snaqplot
net1
writenewick(net1)  # writes to screen, full precision for branch lengths and γ
writenewick(net1, round=true, digits=2)
writenewick(net1, di=true) # γ omitted: for dendroscope
writenewick(net1, "bestnet_h1.tre") # writes to file: creates or overwrites file
rm("bestnet_h1.tre") # hide
```

From a set of candidate networks, one might simply need to score of each network
to pick the best. Here, the score is the negative log pseudo-deviance, and the
lower the better. See the section to get the score of [Candidate networks](@ref).

### Choosing the number of hybridizations (`hmax`)

We change the `hmax` argument to change the search space to let the network have up to 2 or 3 hybrid nodes:
```julia
net2 = snaq!(net1,raxmlCF, hmax=2, filename="net2", seed=3456)
net3 = snaq!(net0,raxmlCF, hmax=3, filename="net3", seed=4567)
```
and plot them (the optimized networks are identical and they both have a single reticulation):
```@example snaqplot
R"svg(name('snaqplot_net23.svg'), width=7, height=3)" # hide
using RCall                  # to be able to tweak our plot within R
R"layout(matrix(1:2, 1, 2))" # to get 2 plots into a single figure: 1 row, 2 columns
R"par"(mar=[0,0,1,0])        # for smaller margins
plot(net2, showgamma=true);
R"mtext"("hmax=2")           # add text annotation: title here
plot(net3, showgamma=true);
R"mtext"("hmax=3")
R"dev.off()"; # hide
nothing # hide
```
![net23](../assets/figures/snaqplot_net23.svg)

with this screen output for net2 (only 1 hybrid node found):

    MaxNet is (C,D,((B,(A)#H7:::0.804411606649347):10.0,(O,(#H7:::0.19558839335065303,E):0.3135243143217013):0.664066456871298):10.0);
    with -loglik 28.31506721890957

and this output for net3 (again, only 1 hybrid found):

    MaxNet is (D,C,((O,(E,#H7:::0.19558839257941849):0.3135243301652981):0.6640664138384673,(B,(A)#H7:::0.8044116074205815):10.0):10.0);
    with -loglik 28.315067218909626

Each network has a `loglik` attribute, which is its log pseudo deviance:
a multiple of the negative log-likelihood up to a constant (the constant is
such that the score is 0 if the network fits the data perfectly).
The lower the better. We can plot these scores across hybrid values:
```@example snaqplot
scores = [loglik(net0), loglik(net1), loglik(net2), loglik(net3)]
hmax = collect(0:3)
R"svg(name('snaqplot_scores_heuristic.svg'), width=4, height=3)" # hide
R"par"(mar=[2.5,2.5,.5,.5], mgp=[1.4,.4,0], tck=-0.02, las=1, lab=[3,5,7]);  # hide
R"plot"(hmax, scores, type="b", ylab="network score", xlab="hmax", col="blue");
R"dev.off()"; # hide
nothing # hide
```
![scores_heuristic](../assets/figures/snaqplot_scores_heuristic.svg)

Here the slope heuristic suggests a single hybrid node:
the score does not get much better beyond h=1.

We made the plot via R above. A more Julian way would use a Julia plotting
package such as [Gadfly](http://gadflyjl.org/stable/) or
[Plots](http://docs.juliaplots.org/latest/ecosystem/), like this for instance:
```julia
using Gadfly
plot(x=collect(0:3), y=scores, Geom.point, Geom.line)
```

Note that since SNaQ assumes level-1 networks (i.e., no intersecting cycles), it might not be possible to add more hybridizations
to networks with few taxa: 

![level1](../assets/level1.png)


Further, a cool [blog](http://avt.im/blog/2018/03/23/R-packages-ggplot-in-julia) about using ggplot within julia.

### Suggestions and best practices


- For big datasets, do only one run (i.e., setting `runs=1`) first to get an idea of computing time:
  ```julia
  net1 = snaq!(net0,raxmlCF, hmax=3, filename="net3", runs=1)
  ```
- Increase the number of hybridizations sequentially:
  `hmax=0,1,2,...`, and use the best network at `h-1` as starting
  point to estimate the best network at `h`.
- Whenever possible, do many independent runs (default `runs=10`)
  to avoid getting stuck on local maxima.
  We do not need to do all runs sequentially.
  We can parallize by doing different runs on different cores (or computers).
  If we are on a machine or on a cluster that has many different cores,
  we can ask Julia to use these multiple cores, and `snaq!` will send
  different runs to different cores, like we did earlier.
- For long jobs, run as a script in the terminal: `julia runSNaQ.jl`,
  arguments to the script are passed to Julia as a vector called `ARGS`.
  See the example script `runSNaQ.jl` in the folder `data_results/scripts/`.
  more on this topic in here: [Parallel computations](@ref)
