# Candidate networks

## Optimizing parameters for a given network

For a given network topology, we can optimize the branch lengths and
inheritance probabilities (γ) with the pseudolikelihood.
This is useful if we have a few candidate networks to compare.
Each network can be optimized individually, and the network with the best
pseudolikelihood can be chosen.

The score being optimized is a pseudo-deviance, i.e.
a multiple of the negative log pseudo-likelihood up to an additive constant
(the lower the better; a pseudo-deviance of 0 corresponds to a perfect fit).

Following our example in [Estimating a network](@ref),
we can optimize parameters on the true network
(the one originally used to simulate the data). 
Given a table of CFs and a network,
the function [`topologymaxQpseudolik!`](@ref)
returns the same network topology
but with optimized branch lengths and inheritance values:

```@setup fixednetworkoptim
using PhyloNetworks, SNaQ, DataFrames
using Logging # to suppress info messages below
baselogger = global_logger()
mkpath("../assets/figures")
exampledir = joinpath(dirname(pathof(SNaQ)), "..","examples")
raxmltrees = joinpath(exampledir,"raxmltrees.tre")
raxmlCF = readtableCF(DataFrame(tablequartetCF(countquartetsintrees(readmultinewick(raxmltrees), showprogressbar=false)...)))
```

```@repl fixednetworkoptim
truenet = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
net1alt = topologymaxQpseudolik!(truenet, raxmlCF);
writenewick(net1alt, round=true)
loglik(net1alt) # pseudo deviance actually: the lower the better
```
```@example fixednetworkoptim
using PhyloPlots, RCall
R"name <- function(x) file.path('..', 'assets', 'figures', x)" # hide
R"svg(name('truenet_opt.svg'), width=4, height=4)" # hide
R"par"(mar=[0,0,0,0])
plot(net1alt, showgamma=true);
R"dev.off()" # hide
nothing # hide
```
![truenet_opt](../assets/figures/truenet_opt.svg)

We get a score of 29.786,
which is comparable to the score of the SNaQ network (net1: 28.315),
especially compared to the score of the best tree (net0: 53.532).
This begs the question: is the true network within the "range" of uncertainty?
We can run a [Bootstrap](@ref) analysis to measure uncertainty
in our network inference.

For a more thorough optimization, we could change the arguments for the tolerances
(`ftolRel` and `xtolAbs`) used to determine when the parameters are optimized and
the search stops (but the optimization will take longer).
It makes no difference on this small data set.
```julia
net1par = topologymaxQpseudolik!(truenet, raxmlCF, ftolRel=1e-10, xtolAbs=1e-10)
loglik(net1par) # pseudo deviance, actually: the lower the better
```

## Network Score with no optimization

For a network with given branch lengths and γ heritabilies,
we can compute the pseudolikelihood (well, a pseudo-deviance) with:
```@repl fixednetworkoptim
topologyQpseudolik!(truenet,raxmlCF);
loglik(truenet) # again, pseudo deviance
```
This function is not maximizing the pseudolikelihood, it is simply computing the
pseudolikelihood (or deviance) for the given branch lengths and probabilities of
inheritance. At the moment, both of these functions require that the
given network is of level 1 (cycles don't overlap).

## Candidate networks compatible with a known outgroup

If the network was estimated via [`snaq!`](@ref), it might turn out to be impossible
to root our estimated network with a known outgroup.
At this time, `snaq!` does not impose any rooting constraint on the network:
the search for the lowest score considers all level-1 networks, including those
that are incompatible with a known outgroup.
(The monophyly of outgroups is not imposed either, like in many other methods.)

If the estimated network cannot be rooted with the known outgroup,
we can check the `.networks` output file for a possible alternative network.
It has a list of networks that are slight modifications of the best network,
where the modifications changed the direction of one reticulation at a time.
For each modified network, the score was calculated. So if we find in this list
a modified network that has a score close to that of the best network,
and that can be re-rooted with our known root position, then this modified network
is a better candidate than the network with the best score.

Below is what the `net1.networks` file looks like, after performing
the analysis in the section [Network estimation](@ref).
Scroll to the right to see the scores.

    (C,D,((O,(E,#H7:::0.19558838614943078):0.31352437658618976):0.6640664399202987,(B,(A)#H7:::0.8044116138505693):10.0):10.0);, with -loglik 28.31506721890958 (best network found, remaining sorted by log-pseudolik; the smaller, the better)
    (C,D,((O,(E)#H7:::0.8150784689693145):0.9336405757682176,(B,(A,#H7:::0.18492153103068557):0.25386142779877724):1.8758156446611114):10.0);, with -loglik 31.535560380783814
    (B,#H7:9.90999345612101::0.2555404440833535,(A,(E,(O,((C,D):10.0)#H7:0.3419231810962026::0.7444595559166465):0.19994859441332047):2.5014911511063644):0.7957621793330066);, with -loglik 56.64548310161462
    (C,D,((O,(E,((B)#H7:::0.7957543284159452,A):4.786202415937916):0.004527712280136759):1.7952610454570868,#H7:::0.20424567158405482):10.0);, with -loglik 67.17775727492258
    (C,D,(#H7:::0.32947301811471164,(B,(A,(E,(O)#H7:::0.6705269818852884):1.371799259141243):0.0):6.397073999864152):7.677245926003807);, with -loglik 199.11401961057143

We can read this file and look at its list of networks like this:

```@repl fixednetworkoptim
file = "net1.networks";
# or use the example file available with the package:
file = joinpath(dirname(pathof(SNaQ)), "..","examples","net1.networks");
netlist = readmultinewick(file) # read the full list of networks in that file
```
Next, we would like to extract the network scores from the file.
Below is some code for doing this in julia
```@repl fixednetworkoptim
score_in_string = read(file, String); # read the file as a single string
score_in_string = (x-> x[1]).(collect(eachmatch(r"with -loglik ([0-9]+.[0-9]+)",score_in_string))); # find all occurences of the loglik scores
scores = parse.(Float64, score_in_string); # parse those matches into numbers

# next: update the "loglik" of each network with the score read from the file
for i in eachindex(netlist)
   netlist[i].fscore = scores[i]
   println("net $i in the list: score = ",scores[i])
end
```
The first network in the list is the best network returned by [`snaq!`](@ref).
We see that the second network has a score that's not too far, but the other networks
have worse scores. The best network and its best modification (second network in the
list) are shown below. We chose to show edge numbers, to use them later
to re-root the networks.

```@example fixednetworkoptim
R"svg(name('fixednetworkoptim_othernets1.svg'), width=7, height=4)" # hide
R"layout(matrix(1:2,1,2))"; # hide
R"par"(mar=[0,0,0,0]) # hide
plot(netlist[1], showgamma=true, showedgenumber=true, tipoffset=0.1);
R"mtext"("best net, score=28.3", line=-1);
plot(netlist[2], showgamma=true, showedgenumber=true, tipoffset=0.1);
R"mtext"("direction modified, score=31.5", line=-1);
R"dev.off()"; # hide
nothing # hide
```
![othernets before reroot](../assets/figures/fixednetworkoptim_othernets1.svg)

Now imagine that our outgroup is taxon A.

Best network: we would get a `RootMismatch` error if we tried to set
the root on the external edge 9 to A, with `rootatnode!(netlist[1], "A")`
(see the PhyloNetworks guide
[Does the root conflict with the direction of a reticulation?](@extref)).
But we could root the best network on the major parent edge to A, edge 10
(rooted network on the left below).

```@example fixednetworkoptim
R"svg(name('fixednetworkoptim_othernets2.svg'), width=7, height=7)" # hide
R"layout(matrix(c(1,4,2,3),2,2))"; # hide
R"par"(mar=[0,0,0.5,0]) # hide
rootonedge!(netlist[1], 10); # root best net to make A outgroup
rotate!(netlist[1], -4); # to 'un-cross' edges
rotate!(netlist[1], -6);
rotate!(netlist[1], -5);
plot(netlist[1], showgamma=true, tipoffset=0.1);
R"mtext"("best net, score=28.3", line=-1);
global_logger(NullLogger()); # hide
rootatnode!(netlist[2], "A"); # net with modified direction: first way to make A outgroup
global_logger(baselogger);   # hide
rotate!(netlist[2], -4) # to 'un-cross' edges
rotate!(netlist[2], -6)
plot(netlist[2], showgamma=true, tipoffset=0.1);
R"mtext"("second best in list, score=31.5\nrequires unsampled population", line=-2);
rootonedge!(netlist[2], 10) # net with modified direction: second way to make A outgroup
for i in [9,-7] rotate!(netlist[2], i); end; # to 'un-cross' edges
plot(netlist[2], showgamma=true, tipoffset=0.1);
R"mtext"("second best in list, score=31.5\ndifferent root position", line=-2);
R"dev.off()"; # hide
nothing # hide
```
![othernets after reroot](../assets/figures/fixednetworkoptim_othernets2.svg)

For the second best network in our list, there are 2 ways to root it with A.
These 2 options give quite different rooted versions of the network:
1. On the external edge 8 to A (top right).
   This requires the existence of an unsampled taxon,
   sister to BDCOE, that would have contributed to introgression into
   an ancestor of E.
2. On its parent edge 10 (bottom right).

A is an outgroup in both rootings, but the second option is more parsimonious,
in the sense that it does not outright require the existence of a "ghost"
taxon: a taxon that went extinct after the introgression, or that is unsampled.

This second rooted version is consistent with 2 possibilities.
It could arise either from  
(a) an unsampled taxon sister to A that contributed to introgression, or  
(b) a direct ancestor of A could have contributed to the
introgression into the ancestor of E.  
Case (a) stipulates the existence of an unsampled ("ghost") taxon,
but case (b) does not require any unsampled taxon.

These two possibilities differ in the length of their gene flow edge
(light blue, here with γ=0.185).
If gene flow came from an unsampled taxon under case (a),
this edge would have positive length, in calendar time.
If gene flow came from a direct ancestor of A under case (b),
the gene flow edge would have length 0.

Edge lengths from SNaQ should be interpreted with caution to
distinguish between the two possibilities because:
* edge lengths estimated with SNaQ are in coalescent units instead of
  calendar time, and necessarily include estimation error;
* an edge with a true length of 0 may be estimated to have a non-zero length
  in coalescent units due to errors in estimated gene trees, to help explain
  gene tree discordance;
* an incorrect topology may result in edges of estimated length 0;
* and some edge lengths are not identifiable from quartet concordance
  factors anyway.

To distinguish between these possibilities, models that separate calendar time,
substitution rate, and population size can be useful. Using stronger assumptions
than SNaQ (e.g. a molecular clock), the rooted network and branch lengths may
be identifiable, so networks with / without unsampled taxa may be distinguished.
See for example [Zhang et al. 2024](https://doi.org/10.1111/tpj.16859),
who used BPP. Note that the model named
"[MSci](https://bpp.github.io/bpp-manual/bpp-4-manual/#the-msc-i-model)"
in BPP is exactly the same as the network multispecies coalescent,
and is typically named NMSC in most papers.
