```@setup expCFs
using PhyloNetworks, SNaQ, DataFrames
mkpath("../assets/figures")
exampledir = joinpath(dirname(pathof(SNaQ)), "..","examples")
raxmltrees = joinpath(exampledir,"raxmltrees.tre")
raxmlCF = readtableCF(DataFrame(tablequartetCF(countquartetsintrees(readmultinewick(raxmltrees), showprogressbar=false)...)))
truenet = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
```

# Quartet test for goodness of fit

One can formally assess whether an estimated network fits the concordance factor data
with the [`QuartetNetworkGoodnessFit`](https://github.com/JuliaPhylo/QuartetNetworkGoodnessFit.jl) package in Julia.
This package simulates concordance factors from the estimated network under the
network multispecies coalescent,
computes how often outlier quartet topologies are observed,
and compares this to the observed concordance factors to
perform a test for goodness-of-fit.
See the [`QuartetNetworkGoodnessFit` documentation](https://juliaphylo.github.io/QuartetNetworkGoodnessFit.jl/stable/man/gof/)
for an example of this test.


# Visualizing observed and expected CFs

A good way to visualize the "goodness-of-fit" of a given estimated network to the data
is to plot the observed CF versus the expected CF. If the network is a good fit, then the dots
in the plot will be close to the diagonal (x=y line).

The following function will create a dataframe with the observed and expected CFs,
which are all saved in the DataCF object after running snaq:
```@repl expCFs
topologymaxQpseudolik!(truenet, raxmlCF);
df_wide = fittedquartetCF(raxmlCF) # same as fittedquartetCF(raxmlCF, :wide)
df_long = fittedquartetCF(raxmlCF, :long)
```
It is important to have run `snaq!`, `topologyQpseudolik!` or `topologymaxQpseudolik!`
before making these tables, or the result would be meaningless.
These functions update the fitted concordance factors (those expected under the network)
inside the DataCF object `raxmlCF`.

Here is one way to plot them, via R again, and using the R package `ggplot2`.

```@example expCFs
using RCall
obsCF = df_long[!,:obsCF]; expCF = df_long[!,:expCF]; # hide
R"name <- function(x) file.path('..', 'assets', 'figures', x)"; # hide
R"svg(name('expCFs_obsvsfitted.svg'), width=5, height=4)"; # hide
R"par"(mar=[2.5,2.6,.5,.5], mgp=[1.5,.4,0], tck=-0.01, las=1, pty="s"); # hide
R"plot(0:1, 0:1, type='l', bty='L', lwd=0.3, col='#008080', xlab='quartet CF observed in gene trees', ylab='quartet CF expected from network')"; # hide
R"set.seed"(1234); # hide
R"points(jitter($obsCF,amount=0.005),jitter($expCF,amount=0.005),col='#008080',bg='#00808090',pch=21)"; # hide
R"dev.off()"; # hide
nothing # hide
```
To install ggplot2 if not installed already, do:
`R"install.packages('ggplot2', dep=TRUE)"`

```julia
@rlibrary ggplot2
ggplot(df_long, aes(x=:obsCF,y=:expCF)) + theme_classic() +
    geom_segment(x=0,y=0,xend=1,yend=1, color="#008080", size=0.3) + # diagonal line
    geom_point(alpha=0.5, color="#008080", position=position_jitter(width=0.005, height=0.005)) +
    ylab("quartet CF expected from network") + xlab("quartet CF observed in gene trees") + coord_equal(ratio=1);
# if needed, save with:
ggsave("expCFs_obsvsfitted.svg", scale=1, width=6, height=5);
```

![obsvsfitted](../assets/figures/expCFs_obsvsfitted.svg)

Many points are overlapping, so they were "jittered" a little to see them all better.
There are always many points overlapping on the bottom-left corner:
concordance factors of 0.0 for quartet resolutions not observed, and not expected.  
To export the table of quartet CFs and explore the fit of the network with other tools:

```julia
using CSV
CSV.write("fittedCF.csv", df_long)
```

We could highlight quartets that include taxon A, say,
if we suspect that it is an unrecognized hybrid.
Many points are overlapping, like before, so they are again "jittered" a bit.

```@example expCFs
using DataFrames
df_long[!,:has_A] .= "no"; # add a column to our data, to indicate which 4-taxon sets have A or not
for r in eachrow(df_long)
    if "A" ∈ [r[:tx1], r[:tx2], r[:tx3], r[:tx4]]
       r[:has_A]="yes"
    end
end
has_A = df_long.has_A # hide
nq = length(has_A); # hide
R"colA=rep('#008080',$nq); bgA=rep('#00808090',$nq);"; # hide
R"colA[$has_A=='yes']='#F8766D'; bgA[$has_A=='yes']='#F8766D90'"; # hide
R"svg(name('expCFs_obsvsfitted_A.svg'), width=5, height=4)"; # hide
R"par"(mar=[2.5,2.6,.5,.5], mgp=[1.5,.4,0], tck=-0.01, las=1, pty="s"); # hide
R"plot(0:1, 0:1, type='l', bty='L', lwd=0.3, col='black', xlab='quartet CF observed in gene trees', ylab='quartet CF expected from network')"; # hide
R"set.seed"(2345); # hide
R"points(jitter($obsCF,amount=0.005),jitter($expCF,amount=0.005),col=colA,bg=bgA,pch=21)"; # hide
R"legend(x=0.7,y=0.3,pch=21,col=c('#008080','#F8766D'),legend=c('no','yes'),title='has A?', bty='n',bg=c('#00808090','#F8766D90'))"; # hide
R"dev.off()"; # hide
nothing # hide
```
```@repl expCFs
first(df_long, 7) # first 7 rows
```

```julia
ggplot(df_long, aes(x=:obsCF, y=:expCF, color=:has_A)) + theme_classic() +
    geom_segment(x=0,y=0,xend=1,yend=1, color="black", size=0.3) + # diagonal line
    geom_point(alpha=0.5, position=position_jitter(width=0.005, height=0.005)) +
    ylab("quartet CF expected from network") + xlab("quartet CF observed in gene trees") + coord_equal(ratio=1);
# can be saved:
ggsave("expCFs_obsvsfitted_A.svg", width=6, height=5);
```

![obsvsfitted A present or not](../assets/figures/expCFs_obsvsfitted_A.svg)
