```@setup multialleles
using PhyloNetworks, SNaQ
```

# Multiple alleles per species

## Between-species 4-taxon sets

The default setting for SNaQ considers that each allele in a gene tree corresponds
to a taxon (a tip) in the network. If instead each allele/individual can be mapped confidently
to a species, and if only the species-level network needs to be estimated,
then the following functions can be used:

```@repl multialleles
using CSV, DataFrames
mappingfile = joinpath(dirname(pathof(SNaQ)), "..","examples",
    "mappingIndividuals.csv");
tm = CSV.read(mappingfile, DataFrame) # taxon map as a data frame
taxonmap = Dict(r[:individual] => r[:species] for r in eachrow(tm)) # as dictionary
```

The [mapping file](https://github.com/juliaphylo/SNaQ/blob/main/examples/mappingIndividuals.csv)
can be a text (or `csv`) file with two columns (at least):
one for the individuals, named `allele` or `individual`,
and one column containing the species names, named `species`.
Each row should map an allele name to a species name.
Next, read in the [gene trees](https://github.com/juliaphylo/SNaQ/blob/main/examples/genetrees_alleletips.tre)
and calculate the quartet CFs at the species level:


```@repl multialleles
genetreefile = joinpath(dirname(pathof(SNaQ)), "..","examples",
    "genetrees_alleletips.tre");
genetrees = readmultinewick(genetreefile);
sort(tiplabels(genetrees[1])) # multiple tips in species S1
df_sp = tablequartetCF(countquartetsintrees(genetrees, taxonmap;
    showprogressbar=false)...);
keys(df_sp)  # columns names
df_sp[:qind] # quartet index
df_sp[:t1]   # name of first taxon in each quartet
df_sp[:CF12_34] # concordance factor for split taxa 12 vs 34
```

Now `df_sp` is a table (technically, a named tuple) containing the
quartet concordance factors at the species level only, that is,
considering sets made of 4 distinct species,
even if the gene trees may have multiple alleles from the same species.
For 4 distinct species `A,B,C,D`, all alleles from each species (`A` etc.)
will be used to calculate the quartet CF. If a given gene tree has
`n_a` alleles from `a`, `n_b` alleles from `b` etc., then
each set of 4 alleles is given a weight of `1/(n_a n_b n_c n_d)`
to calculated of the CF for `A,B,C,D` (such that the total weight from
this particular gene trees is 1).
It is safe to save this data frame, then use it for `snaq!` like this:

```@repl multialleles
CSV.write("tableCF_species.csv", df_sp);   # to save the table to a file
d_sp = readtableCF("tableCF_species.csv"); # "DataCF" object for use in snaq!
summarizedataCF(d_sp)
```

## Within-species 4-taxon sets

Four-taxon sets involving 2 individuals per species can provide more
information about the underlying network, including external branch
length in coalescent units. However, [`snaq!`](@ref) runs more slowly when
using this extra information.

To get quartet CFs from sets of 4 individuals
in which 2 individuals are from the same species, the following functions
should be used,
where the mapping file can be a text (or `csv`) file with two columns
named `allele` (or `individual`) and `species`,
mapping each allele name to a species name.

```@repl multialleles
q,t = countquartetsintrees(genetrees);
df_ind = DataFrame(tablequartetCF(q,t)); # no mapping: CFs across individuals
first(df_ind, 5) # to see the first 5 rows
```
Now `df_ind` is the table of concordance factors at the level of individuals.
In other words, it lists CFs using one row for each set of 4 alleles/individuals.

**Warning**:
This procedure requires that all alleles from the same individual are given
the same name (the individual's 'name') *across all genes* for which that
individual was sequenced.

Next, we use [`mapallelesCFtable`](@ref) to get these data as
quartet concordance factors at the species level in `df_sp`:
with the allele names replaced by the appropriate species names.

```@repl multialleles
CSV.write("tableCF_individuals.csv", df_ind);  # to save to a file
df_sp = mapallelesCFtable(mappingfile, "tableCF_individuals.csv";
    columns=2:5); # taxon names are in columns 2 through 5, not default 1-4
nrow(df_sp)       # 35 quartets of individuals
first(df_sp, 6)   # first 6 rows of data frame
```

The warning above, after creating `df_sp`, is because our mapping file does not
list species S2 through S5. We did not need to list them because we have a
single individual in each of these species.
So we can safely ignore the warning.
We will just need to make sure that our starting tree, when we run SNaQ,
has the same (unmapped) names, here S2-S5.

The command below modifies `df_sp` to delete rows that are uninformative about
between-species relationships, such as rows containing 3 or 4 individuals from
the same species (e.g. rows 1, 2 and 6: they contain S1 three times);
and creates `d_sp`, an object of type [`DataCF`](@ref) at the species level,
that we can use later as input for networks estimation with [`snaq!`](@ref).
```@repl multialleles
d_sp = readtableCF!(df_sp, mergerows=true); # DataCF object
nrow(df_sp) # 31 quartets of individuals informative between species
```

!!! note
    For a four-taxon set `A,B,C,D`, all the individuals from `A`, `B`, `C` and `D`
    are considered, say `(a1,b1,c1,d1)`, `(a2,b1,c1,d1)`, `(a1,b2,c1,d1)`,
    `(a2,b2,c1,d1)` and so on.
    The CFs of these 4-taxon sets are averaged together to obtain the
    CFs at the species level. This procedures gives more weight to genes that have
    many alleles (because they contribute to more sets of 4 individuals) and less
    weight to genes that have few alleles.

Before we run SNaQ, it is safe to save the concordance factor of *species* quartets,
which can be calculated by averaging the CFs of quartets of individuals
from the associated species:

```@repl multialleles
df_sp_ave = DataFrame(tablequartetCF(d_sp))  # CFs averaged across individuals
CSV.write("CFtable_species.csv", df_sp_ave); # save to file
```

Some quartets have the same species repeated twice,
representing cases when 2 of the 4 individuals came from the same species.
These quartets, with repeated species, are informative about the population
size of extant populations, i.e. about the lengths of external branches in
coalescent units.

The main difference between this section compared to the previous section
on [Between-species 4-taxon sets](@ref) is that quartets with 2 individuals from
the same species are included here, such as `a1,a2,b1,c1`.
Also, the weighting of quartets is different. Here, genes with more alleles
are given more weight.

now we can run `snaq!`:

```julia
net = snaq!(T_sp, d_sp);
```
where `T_sp` should be a starting topology with one tip per species,
labelled with the same species names as the names used in the mapping file.

If `snaq!` takes too long that way, we can try a less ambitious estimation
that does not estimate the external branch lengths, that is,
*without* using quartets that have 2 individuals from the same species.
To do so, we can use the quartet concordance factors at the species level,
but filter out the quartets with one (or more) species repeated, such as
these below:

```@repl multialleles
first(df_sp_ave, 3) # some quartets have the same species twice
```

Filtering them out can be done as in the first section
([Between-species 4-taxon sets](@ref)) to give equal weight to all genes,
or as shown below to give more weight to genes that have more alleles.
We first define a helper function to identify which rows we want to get rid of.

```@repl multialleles
first(df_sp_ave, 3) # some quartets have the same species twice
"""
    hasrep

Return true if a row (4-taxon set) has a "repeated" species, that is, a species
whose name ends with "__2". Otherwise, return false.

Warning: this function assumes that taxon names are in columns
"t1", "t2", "t3", "t4". For data frames with different column names,
e.g. "taxon1", "taxon2" etc., simply edit the code below by replacing
`:t1` by `:taxon1` (or the appropriate column name in your data).
"""
function hasrep(row)
    occursin(r"__2$", row[:t1]) || occursin(r"__2$", row[:t2]) ||
    occursin(r"__2$", row[:t3]) || occursin(r"__2$", row[:t4])
end; # this function is used on the next line
df_sp_reduced = filter(!hasrep, df_sp_ave) # removes rows with repeated species
CSV.write("CFtable_species_norep.csv", df_sp_reduced); # to save to file
d_sp_reduced = readtableCF(df_sp_reduced) # DataCF object, for input to snaq!
```

and now we can run `snaq!` on the reduced set of quartets without repeats,
which should be faster:

```julia
net = snaq!(T_sp, d_sp_reduced);
```
