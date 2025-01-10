# function to write a csv table from the expCF of an
# array of quartets
# warning: does not check if the expCF have been calculated
function writeExpCF(quartets::Array{Quartet,1})
    df = DataFrames.DataFrame(t1=String[],t2=String[],t3=String[],t4=String[],
                              CF12_34=Float64[],CF13_24=Float64[],CF14_23=Float64[])
    for q in quartets
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.qnet.expCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        push!(df, [q.taxon[1],q.taxon[2],q.taxon[3],q.taxon[4],q.qnet.expCF[1],q.qnet.expCF[2],q.qnet.expCF[3]])
    end
    return df
end
writeExpCF(d::DataCF) = writeExpCF(d.quartet)

"""
    tablequartetCF(vector of Quartet objects)
    tablequartetCF(DataCF)
    tablequartetCF(gene trees)
    tablequartetCF(gene tree file)

Build a NamedTuple containing observed quartet concordance factors,
with the fields named:
- `t1`, `t2`, `t3`, `t4` for the four taxon names in each quartet
- `CF12_34`, `CF13_24`, `CF14_23` for the 3 quartets of a given four-taxon set
- `ngenes` if this information is available for some quartets

Some downstream functions may require the observed quartet concordance factors
to be in a DataFrame, this can be easily converted by wrapping the output
NamedTuple in the `DataFrame()` function
"""
function tablequartetCF(quartets::Array{Quartet,1})
    nq = length(quartets)

    nt  = (t1=Array{String}(undef,nq),
            t2=Array{String}(undef,nq),
            t3=Array{String}(undef,nq),
            t4=Array{String}(undef,nq),
            CF12_34=Array{Float64}(undef,nq),
            CF13_24=Array{Float64}(undef,nq),
            CF14_23=Array{Float64}(undef,nq),
            ngenes=Array{Union{Missing, Float64}}(undef,nq))

    for (i,q) in enumerate(quartets)
        length(q.taxon) == 4 || error("quartet $(q.number) does not have 4 taxa")
        length(q.obsCF) == 3 || error("quartet $(q.number) does have qnet with 3 expCF")
        
        nt.t1[i] = q.taxon[1]
        nt.t2[i] = q.taxon[2]
        nt.t3[i] = q.taxon[3]
        nt.t4[i] = q.taxon[4]
        nt.CF12_34[i] = q.obsCF[1]
        nt.CF13_24[i] = q.obsCF[2]
        nt.CF14_23[i] = q.obsCF[3]
        nt.ngenes[i] = (q.ngenes==-1.0 ? missing : q.ngenes)
        
    end
    if all(ismissing, nt.ngenes)
        nt = Base.structdiff(nt,NamedTuple{(:ngenes,)}((nothing,)))
    end
    return nt
end
tablequartetCF(d::DataCF) = tablequartetCF(d.quartet)




"""
    readtableCF(file)
    readtableCF(data frame)
    readtableCF!(data frame)

Read a file or DataFrame object containing a table of concordance factors (CF),
with one row per 4-taxon set. The first 4 columns are assumed to give the labels
of the 4 taxa in each set (tx1, tx2, tx3, tx4).
Columns containing the CFs are assumed to be named
`CF12_34`, `CF13_24` and `CF14_23`;
or `CF12.34`, `CF13.24` and `CF14.23`;
or else are assumed to be columns 5,6,7.
If present, a column named 'ngenes' will be used to get the number of loci
used to estimate the CFs for each 4-taxon set.

Output: [`DataCF`](@ref) object

Optional arguments:

- `summaryfile`: if specified, a summary file will be created with that name.
- `delim` (for the first form only): to specify how columns are delimited,
  with single quotes: delim=';'. Default is a `csv` file, i.e. `delim=','`.
- `mergerows`: false by default. When true, will attempt to merge multiple rows
  corresponding to the same four-taxon set (by averaging their quartet CFs) even
  if none of the species is repeated within any row (that is, in any set of 4 taxa)

The last version modifies the input data frame, if species are represented by multiple alleles
for instance (see [`readtableCF!`](@ref)(data frame, columns)).
"""
function readtableCF(file::AbstractString; delim=','::Char, summaryfile=""::AbstractString, kwargs...)
    df = DataFrame(CSV.File(file, delim=delim); copycols=false)
    readtableCF!(df; summaryfile=summaryfile, kwargs...)
end

function readtableCF(df0::DataFrames.DataFrame; summaryfile=""::AbstractString, kwargs...)
    df = deepcopy(df0)
    readtableCF!(df; summaryfile=summaryfile, kwargs...)
end

function readtableCF!(df::DataFrames.DataFrame; summaryfile=""::AbstractString, kwargs...)
    @debug "assume the numbers for the taxon read from the observed CF table match the numbers given to the taxon when creating the object network"
    taxoncolnames = [[:t1, :tx1, :tax1, :taxon1], [:t2, :tx2, :tax2, :taxon2],
                     [:t3, :tx3, :tax3, :taxon3], [:t4, :tx4, :tax4, :taxon4] ]
    taxoncol = [findfirst(x-> x ∈ taxoncolnames[1], DataFrames.propertynames(df)),
                findfirst(x-> x ∈ taxoncolnames[2], DataFrames.propertynames(df)),
                findfirst(x-> x ∈ taxoncolnames[3], DataFrames.propertynames(df)),
                findfirst(x-> x ∈ taxoncolnames[4], DataFrames.propertynames(df))]
    alternativecolnames = [ # obsCF12 is as exported by fittedquartetCF()
        [:CF12_34, Symbol("CF12.34"), :obsCF12, :CF1234],
        [:CF13_24, Symbol("CF13.24"), :obsCF13, :CF1324],
        [:CF14_23, Symbol("CF14.23"), :obsCF14, :CF1423]
    ]
    obsCFcol = [findfirst(x-> x ∈ alternativecolnames[1], DataFrames.propertynames(df)),
                findfirst(x-> x ∈ alternativecolnames[2], DataFrames.propertynames(df)),
                findfirst(x-> x ∈ alternativecolnames[3], DataFrames.propertynames(df))]
    ngenecol =  findfirst(isequal(:ngenes), DataFrames.propertynames(df))
    withngenes = ngenecol !== nothing
    nothing in taxoncol && error("columns for taxon names were not found")
    nothing in obsCFcol && error(
      """Could not identify columns with quartet concordance factors (qCFs).
      Was expecting CF12_34, CF13_24 and CF14_23 for the columns with CF values,
      or CF12.34 or obsCF12, etc.""")
    columns = [taxoncol; obsCFcol]
    if withngenes  push!(columns, ngenecol)  end

    d = readtableCF!(df, columns; kwargs...)

    if withngenes # && d.numTrees == -1
        m1 = minimum([q.ngenes for q in d.quartet])
        m2 = maximum([q.ngenes for q in d.quartet])
        if m1<m2 print("between $m1 and ") end
        println("$m2 gene trees per 4-taxon set")
        # other info printed by show() on a DataCF object: num quartets and num gene trees
    end
    if(summaryfile != "")
        descData(d,summaryfile)
    end
    return d
end


# see docstring below, for readtableCF!
# takes in df and 7 or 8 column numbers (4 labels + 3 CFs + ngenes possibly)
function readtableCF!(df::DataFrames.DataFrame, co::Vector{Int}; mergerows=false)
    withngenes = (length(co)==8) # true if column :ngenes exists, false ow
    repSpecies = cleanAlleleDF!(df,co) # removes uninformative rows from df (not df0)
    # fixit: cleanAlleleDF! is time consuming but many times not needed
    # add option to skip it, if the user knows that each tip appears once only?
    if mergerows || !isempty(repSpecies)
        df = mergeRows(df,co)   # warning: this 'df' is *not* changed externally
        co = collect(eachindex(co)) # 1:7 or 1:8
    end                         # we cannot move to mapallelesCFtable because we need repSpecies in here
    quartets = Quartet[]
    for i in 1:size(df,1)
        push!(quartets,Quartet(i,string(df[i,co[1]]),string(df[i,co[2]]),string(df[i,co[3]]),string(df[i,co[4]]),
                               [df[i,co[5]],df[i,co[6]],df[i,co[7]]]))
        if withngenes
            quartets[end].ngenes = df[i,co[8]]
        end
    end
    d = DataCF(quartets)
    if(!isempty(repSpecies))
        d.repSpecies = repSpecies
    end
    return d  # return d, df ## to save memory & gc with readtableCF! for bootstrapping?
end


"""
    readtableCF!(data frame, columns; mergerows=false)

Read in quartet CFs from data frame, assuming information is in columns numbered `columns`,
of length **7 or 8**: 4 taxon labels then 3 CFs then ngenes possibly.

If some species appears more than once in the same 4-taxon set (e.g. t1,t1,t2,t3),
then the data frame is modified to remove rows (4-taxon sets) that are uninformative about
between-species relationships. This situation may occur if multiple individuals are sampled from
the same species. A 4-taxon set is uninformative (and its row is removed)
if one taxon is repeated 3 or 4 times (like t1,t1,t1,t1 or t1,t2,t2,t2).
The list of species appearing twice in some 4-taxon sets is stored in the output DataCF object.
For these species, the length of their external edge is identifiable (in coalescent units).
If multiple rows correspond to the same 4-taxon set, these rows are merged and their CF values
(and number of genes) are averaged.
If none of the species is repeated within any 4-taxon set, then this averaging
is attempted only if `mergerows` is true.

    readtableCF!(DataCF, data frame, columns)

Modify the `.quartet.obsCF` values in the `DataCF` object with those read from the data frame
in columns numbered `columns`.
`columns` should have **3** columns numbers for the 3 CFs in this order:
`12_34`, `13_24` and `14_23`.

Assumptions:
- same 4-taxon sets in `DataCF` and in the data frame, and in the same order,
  but this assumption is *not checked* (for speed, e.g. during bootstrapping).
- one single row per 4-taxon set (multiple individuals representatives
  of the same 4-taxon set should have been already merged);
  basically: the DataCF should have been created from the data frame by `readtableCF!(df, colums)`
"""
function readtableCF!(datcf::DataCF, df::DataFrame, cols::Vector{Int})
    for i in 1:size(df,1)
        for j in 1:3
            datcf.quartet[i].obsCF[j] = df[i,cols[j]]
        end
    end
end


# ---------------- read input gene trees and calculate obsCF ----------------------

"""
    readmultinewicklevel1(file)

Read a text file with a list of trees/networks in extended newick format
(one tree per line) and transform them like [`readnewicklevel1`](@ref).
Namely, in each tree/network
- the root is suppressed (becomes of degree 3 if it was of degree 2)
- any polytomy is resolved arbitrarily
- any missing branch length is set to 1
- any branch length above 10 is set to 10 (this assumes branch lengths in coalescent units)
- any missing γ's are set to (0.1, 0.9)
and more (see [`readnewicklevel1`](@ref)).

See [`PhyloNetworks.readmultinewick`]()
to read multiple trees or networks with no modification.

Output: array of HybridNetwork objects.

Each line starting with "(" will be considered as describing one topology.
The file can have extra lines that are ignored.
"""
function readmultinewicklevel1(file::AbstractString)
    try
        s = open(file)
    catch
        error("Could not find or open $(file) file");
    end
    vnet = HybridNetwork[];
    s = open(file)
    numl = 1
    for line in eachline(s)
        line = strip(line) # remove spaces
        @debug "$(line)"
        c = isempty(line) ? "" : line[1]
        if(c == '(')
           try
               push!(vnet, readTopologyUpdate(line,false))
           catch(err)
               error("could not read tree in line $(numl). The error is $(err)")
           end
        end
        numl += 1
    end
    close(s)
    return vnet # consistent output type: HybridNetwork vector. might be of length 0.
end


# function to list all quartets for a set of taxa names
# return a vector of quartet objects, and if writeFile=true, writes a file
function allQuartets(taxon::Union{Vector{<:AbstractString},Vector{Int}}, writeFile::Bool)
    vquartet = Quartet[];
    if(writeFile)
        #allName = "allQuartets$(string(integer(time()/1000))).txt"
        allName = "allQuartets.txt"
        f = open(allName,"w")
    end
    i = 1
    
    # Iterate over all quartet combinations (this is equivalent to enumerating `combinations(taxon, 4)`
    # but does not require the `Combinatorics` pkg dependency)
    for taxa_idx1 = 1:(length(taxon) - 3)
        for taxa_idx2 = (taxa_idx1+1):(length(taxon)-2)
            for taxa_idx3 = (taxa_idx2+1):(length(taxon)-1)
                for taxa_idx4 = (taxa_idx3+1):length(taxon)
                    if(writeFile)
                        write(f,"$(taxon[taxa_idx1]),$(taxon[taxa_idx2]),$(taxon[taxa_idx3]),$(taxon[taxa_idx4])\n")
                    end
                    push!(vquartet,Quartet(i,string(taxon[taxa_idx1]),string(taxon[taxa_idx2]),string(taxon[taxa_idx3]),chomp(string(taxon[taxa_idx4])),[1.0,0.0,0.0]))
                    i += 1 # overflow error if # quartets > typemax(Int), i.e. if 121,978+ taxa with Int64, 478+ taxa with Int32
                end
            end
        end
    end
    if(writeFile)
        close(f)
    end
    return vquartet
end
allQuartets(numtaxa::Integer, writeFile::Bool) = allQuartets(1:numtaxa, writeFile)


# function to list num randomly selected quartets for a vector of all quartets
# return a vector of Quartet randomly chosen (and a file if writeFile=true)
function randQuartets(allquartets::Vector{Quartet},num::Integer, writeFile::Bool)
    randquartets = Quartet[]
    n = length(allquartets)
    if num == 0
        num = Integer(floor(0.1*n))
    end
    num <= n || error("you cannot choose a sample of $(num) quartets when there are $(n) in total")
    indx = [ones(Bool,num); zeros(Bool,n-num)]
    indx = indx[sortperm(randn(n))]
    if writeFile
        randName = "rand$(num)Quartets.txt"
        println("list of randomly selected quartets in file $(randName)")
        out = open(randName,"w")
    end
    for i in 1:n
        if indx[i]
            if writeFile
                q = allquartets[i].taxon
                write(out,"$(q[1]),$(q[2]),$(q[3]),$(q[4])\n")
            end
            push!(randquartets,allquartets[i])
        end
    end
    if writeFile
        close(out)
    end
    return randquartets
end


# function that will not use randQuartets(list of quartets,...)
# this function uses whichQuartet to avoid making the list of all quartets
# fixit: writeFile is not used. remove the option?
function randQuartets(taxon::Union{Vector{<:AbstractString},Vector{Int}},num::Integer, writeFile::Bool)
    randquartets = Quartet[]
    n = length(taxon)
    ntotal = binom(n,4)
    num <= ntotal || error("you cannot choose a sample of $(num) quartets when there are $(ntotal) in total")
    # indx = [rep(1,num);rep(0,ntotal-num)] # requires much more memory than necessary:
    # indx = indx[sortperm(randn(ntotal))]  # several arrays of size ntotal !!
    # rq = findall(x -> x==1, indx)
    rq = sample(1:ntotal, num, replace=false, ordered=true)
    randName = "rand$(num)Quartets.txt"
    println("list of randomly selected quartets in file $(randName)")
    out = open(randName,"w")
    i = 1
    for q in rq
        qind = whichQuartet(n,q) # vector of int
        quartet = createQuartet(taxon,qind,i)
        write(out,"$(quartet.taxon[1]), $(quartet.taxon[2]), $(quartet.taxon[3]), $(quartet.taxon[4])\n")
        push!(randquartets,quartet)
        i += 1
    end
    close(out)
    return randquartets
end
randQuartets(numtaxa::Integer,num::Integer, writeFile::Bool) = randQuartets(1:numtaxa,num, writeFile)


# function to read list of quartets from a file
# and create Quartet type objects
function readListQuartets(file::AbstractString)
    try
        f = open(file)
    catch
        error("could not open file $(file)")
    end
    f = open(file)
    quartets = Quartet[];
    i = 1
    for line in eachline(f)
        l = split(line,",")
        length(l) == 4 || error("quartet with $(length(l)) elements, should be 4: $(line)")
        push!(quartets,Quartet(i,string(l[1]),string(l[2]),string(l[3]),chomp(string(l[4])),[1.0,0.0,0.0]))
        i += 1
    end
    return quartets
end


"""
    sameTaxa(Quartet, HybridNetwork)

Return `true` if all taxa in the quartet are represented in the network,
`false` if one or more taxa in the quartet does not appear in the network.

warning: the name can cause confusion. A more appropriate name might be
"in", or "taxain", or "taxonsubset", or etc.
"""
function sameTaxa(q::Quartet, t::HybridNetwork)
    for name in q.taxon
        in(name,t.names) || return false
    end
    return true
end


"""
    taxadiff(Vector{Quartet}, network; multiplealleles=true)
    taxadiff(DataCF, network; multiplealleles=true)

Return 2 vectors:

- taxa in at least 1 of the quartets but not in the network, and
- taxa in the network but in none of the quartets.

When `multiplealleles` is true, the taxon names that end with "__2"
are ignored in the quartets: they are not expected to appear in the
networks that users give as input, or get as output.
"""
function taxadiff(quartets::Vector{Quartet}, t::HybridNetwork;
                  multiplealleles=true::Bool)
    tq = tiplabels(quartets)
    secondallele = occursin.(Ref(r"__2$"), tq)
    for i in length(secondallele):-1:1
        secondallele[i] || continue
        basetax = match(r"(.*)__2$", tq[i]).captures[1]
        # if tq[i] = "mouse__2" for instance, then basetax = "mouse"
        if basetax in tq     # some other taxon is "mouse"
            deleteat!(tq, i) # delete "mouse__2" from tq IF "mouse" is present
        end
    end
    tn = tiplabels(t)
    return (setdiff(tq,tn), setdiff(tn,tq))
end

taxadiff(d::DataCF, t::HybridNetwork; multiplealleles=true::Bool) =
    taxadiff(d.quartet, t; multiplealleles=multiplealleles)





# extract & sort the union of taxa of list of quartets
function tiplabels(quartets::Vector{Quartet})
    taxa = reduce(union, q.taxon for q in quartets)
    return sort_stringasinteger!(taxa)
end
tiplabelsTree(file::AbstractString) = tiplabels(readmultinewicklevel1(file))

tiplabels(d::DataCF) = tiplabels(d.quartet)

"""
    calculateObsCFAll!(DataCF, taxa::Union{Vector{<:AbstractString}, Vector{Int}})

Calculate observed concordance factors:
update the `.quartet[i].obsCF` values of the `DataCF` object based on its .tree vector.

    calculateObsCFAll!(vector of quartets, vector of trees, taxa)

Calculate observed concordance factors:
update the `.obsCF` values of the quartets, based on the trees, and returns a new `DataCF` object
with these updated quartets and trees.

    calculateObsCFAll_noDataCF!(vector of quartets, vector of trees, taxa)

update the `.obsCF` values of the quartets based on the trees, but returns nothing.

Warning: all these functions need input trees (without any reticulations: h=0).

See also: [`PhyloNetworks.countquartetsintrees`](@extref), which uses a faster algorithm,
processing each input tree only once.
`calculateObsCFAll_noDataCF!` processes each input tree `# quartet` times.
"""
function calculateObsCFAll!(dat::DataCF, taxa::Union{Vector{<:AbstractString}, Vector{Int}})
    calculateObsCFAll_noDataCF!(dat.quartet, dat.tree, taxa)
end


function calculateObsCFAll!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{<:AbstractString}, Vector{Int}})
    calculateObsCFAll_noDataCF!(quartets, trees, taxa)
    d = DataCF(quartets,trees)
    return d
end


function calculateObsCFAll_noDataCF!(quartets::Vector{Quartet}, trees::Vector{HybridNetwork}, taxa::Union{Vector{<:AbstractString}, Vector{Int}})
    println("calculating obsCF from $(length(trees)) gene trees and for $(length(quartets)) quartets")
    index = 1
    totalq = length(quartets)
    println("Reading in quartets...")
    r = round(1/totalq, digits=2)
    numq = (r > 0.02 ? totalq : 50)
    print("0+")
    for i in 1:numq
        print("-")
    end
    print("+100%")
    println("  ")
    print("  ")
    for q in quartets
        if round(index/totalq, digits=2) > 0.02
            print("*")
            index = 1
        end
        suma = 0
        sum12 = 0
        sum13 = 0
        sum14 = 0
        for t in trees
            isTree(t) || error("gene tree found in file that is a network $(writenewick(t))")
            if sameTaxa(q,t)
                M = tree2Matrix(t,taxa) #fixit: way to reuse M? length(t.edge) will be different across trees
                res = extractQuartetTree(q,M,taxa)
                @debug "res is $(res)"
                if(res == 1)
                    sum12 += 1
                elseif(res == 2)
                    sum13 += 1
                elseif(res == 3)
                    sum14 += 1
                end
                suma += (res == 0) ? 0 : 1
            end
        end
        q.obsCF = [sum12/suma, sum13/suma, sum14/suma]
        q.ngenes = suma
        index += 1
    end
    println("  ")
    return nothing
end



"""
    readInputData(trees, quartetfile, whichQuartets, numQuartets, writetable, tablename, writeQfile, writesummary)
    readInputData(trees, whichQuartets, numQuartets, taxonlist,   writetable, tablename, writeQfile, writesummary)

Read gene trees and calculate the observed quartet concordance factors (CF),
that is, the proportion of genes (and the number of genes) that display each
quartet for a given list of four-taxon sets.

Input:

- `trees`: name of a file containing a list of input gene trees,
  or vector of trees (`HybridNetwork` objects)

Optional arguments (defaults):

- `quartetfile`: name of a file containing a list of quartets, or more precisely,
  a list of four-taxon sets
- `whichQuartets` (`:all`): which quartets to sample.
  `:all` for all of them, `:rand` for a random sample.
- `numQuartets`: number of quartets in the sample.
  default: total number of quartets if `whichQuartets=:all`
  and 10% of total if `whichQuartets=:rand`
- `taxonlist` (all in the input gene trees):
  If `taxonlist` is used, `whichQuartets` will consist of *all* sets of 4 taxa in the `taxonlist`. 
- `writetable` (true): write the table of observed CF?
- `tablename` ("tableCF.txt"): if `writetable` is true, the table of observed CFs is write to file `tablename`
- `writeQfile` (false): write intermediate file with sampled quartets?
- `writesummary` (true): write a summary file?
  if so, the summary will go in file "summaryTreesQuartets.txt".

Uses [`calculateObsCFAll!`](@ref), which implements a slow algorithm.

See also:
[`PhyloNetworks.countquartetsintrees`](@extref), which uses a much faster algorithm;
[`readtrees2CF`](@ref), which is basically a re-naming of `readInputData`, and
[`readtableCF`](@ref) to read a table of quartet CFs directly.
"""
function readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readtableCF(\"$(filename)\")""")
        end
    end
    println("read input trees from file $(treefile)\nand quartetfile $(quartetfile)")
    trees = readmultinewicklevel1(treefile)
    readInputData(trees, quartetfile, whichQ, numQ, writetab, filename, writeFile, writeSummary)
end
readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool) = readInputData(treefile, quartetfile, whichQ, numQ, writetab, "none", false, true)
readInputData(treefile::AbstractString, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer) = readInputData(treefile, quartetfile, whichQ, numQ, true, "none", false, true)
readInputData(treefile::AbstractString, quartetfile::AbstractString, writetab::Bool, filename::AbstractString) = readInputData(treefile, quartetfile, :all, 0, writetab, filename, false, true)


function readInputData(trees::Vector{HybridNetwork}, quartetfile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || @warn "set numQ=$(numQ) but whichQ is not rand, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option"
        println("will use all quartets in file $(quartetfile)")
        quartets = readListQuartets(quartetfile)
    elseif(whichQ == :rand)
        if(numQ == 0)
            @warn "not specified numQ but whichQ=rand, so 10% of quartets will be sampled" #handled inside randQuartets
        else
            println("will take a random sample of $(numQ) 4-taxon sets from file $(quartetfile)")
        end
        allquartets = readListQuartets(quartetfile)
        quartets = randQuartets(allquartets,numQ,writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees, tiplabels(trees))
    if(writetab)
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readtableCF(\"$(filename)\")""")
        end
        println("\ntable of obsCF printed to file $(filename)")
        df = tablequartetCF(d)
        CSV.write(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    writeSummary && descData(d,"summaryTreesQuartets.txt")
    return d
end


function readInputData(treefile::AbstractString, whichQ::Symbol=:all, numQ::Integer=0,
        taxa::Union{Vector{<:AbstractString}, Vector{Int}}=tiplabelsTree(treefile),
        writetab::Bool=true, filename::AbstractString="none",
        writeFile::Bool=false, writeSummary::Bool=true)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt" # "tableCF$(string(integer(time()/1000))).txt"
        end
        if (isfile(filename) && filesize(filename) > 0)
           error("""file $(filename) already exists and is non-empty. Cannot risk to erase data.
                    Choose a different CFfile name, use writeTab=false, or read the existing file
                    with readtableCF(\"$(filename)\")""")
        end
    end
    println("read input trees from file $(treefile). no quartet file given.")
    trees = readmultinewicklevel1(treefile)
    readInputData(trees, whichQ, numQ, taxa, writetab, filename, writeFile, writeSummary)
end
readInputData(treefile::AbstractString, whichQ::Symbol, numQ::Integer, writetab::Bool) = readInputData(treefile, whichQ, numQ, tiplabelsTree(treefile), writetab, "none",false, true)
readInputData(treefile::AbstractString, taxa::Union{Vector{<:AbstractString}, Vector{Int}}) = readInputData(treefile, :all, 0, taxa, true, "none",false, true)
# above: the use of tiplabelsTree to set the taxon set
#        is not good: need to read the tree file twice: get the taxa, then get the trees
#        this inefficiency was fixed in readtrees2CF


function readInputData(trees::Vector{HybridNetwork}, whichQ::Symbol, numQ::Integer, taxa::Union{Vector{<:AbstractString}, Vector{Int}}, writetab::Bool, filename::AbstractString, writeFile::Bool, writeSummary::Bool)
    if(whichQ == :all)
        numQ == 0 || @warn "set numQ=$(numQ) but whichQ=all, so all quartets will be used and numQ will be ignored. If you want a specific number of 4-taxon subsets not random, you can input with the quartetfile option"
        quartets = allQuartets(taxa,writeFile)
        println("will use all quartets on $(length(taxa)) taxa")
    elseif(whichQ == :rand)
        if(numQ == 0)
            @warn "not specified numQ with whichQ=rand, so 10% of quartets will be sampled" #handled inside randQuartets
        else
            println("will use a random sample of $(numQ) 4-taxon sets ($(round((100*numQ)/binomial(length(taxa),4), digits=2)) percent) on $(length(taxa)) taxa")
        end
        quartets = randQuartets(taxa,numQ, writeFile)
    else
        error("unknown symbol for whichQ $(whichQ), should be either all or rand")
    end
    d = calculateObsCFAll!(quartets,trees,taxa)
    if writetab
        if(filename == "none")
            filename = "tableCF.txt"
        end
        println("table of obsCF printed to file $(filename)")
        df = tablequartetCF(d)
        CSV.write(filename,df)
    end
    #descData(d,"summaryTreesQuartets$(string(integer(time()/1000))).txt")
    writeSummary && descData(d,"summaryTreesQuartets.txt")
    return d
end



# rename the function readInputData to make it more user-friendly
"""
    readtrees2CF(treefile)
    readtrees2CF(vector of trees)

Read trees in parenthetical format from a file, or take a vector of trees already read,
and calculate the proportion of these trees having a given quartet (concordance factor: CF),
for all quartets or for a sample of quartets.
Optional arguments include:

- `quartetfile`: name of text file with list of 4-taxon subsets to be analyzed. If none is specified, the function will list all possible 4-taxon subsets.
- `whichQ="rand"`: to choose a random sample of 4-taxon subsets
- `numQ`: size of random sample (ignored if whichQ is not set to "rand")
- `writeTab=false`: does not write the observedCF to a table (default true)
- `CFfile`: name of file to save the observedCF (default tableCF.txt)
- `writeQ=true`: save intermediate files with the list of all 4-taxon subsets and chosen random sample (default false).
- `writeSummary`: write descriptive stats of input data (default: true)
- `nexus`: if true, it assumes the gene trees are written in nexus file (default: false)

See also:
[`PhyloNetworks.countquartetsintrees`](@extref), which uses a much faster algorithm;
[`readtableCF`](@ref) to read a table of quartet CFs directly.
"""
function readtrees2CF(treefile::AbstractString; quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer,
                      writeTab=true::Bool, CFfile="none"::AbstractString,
                      taxa::AbstractVector=Vector{String}(),
                      writeQ=false::Bool, writeSummary=true::Bool, nexus=false::Bool)
    trees = (nexus ?
             readnexus_treeblock(treefile, readTopologyUpdate, false, false; reticulate=false) :
             readmultinewicklevel1(treefile))
    if length(taxa)==0        # tiplabels(trees) NOT default argument:
      taxa = tiplabels(trees) # otherwise: tree file is read twice
    end
    readtrees2CF(trees, quartetfile=quartetfile, whichQ=whichQ, numQ=numQ, writeTab=writeTab,
                 CFfile=CFfile, taxa=taxa, writeQ=writeQ, writeSummary=writeSummary)
end


# same as before, but with input vector of HybridNetworks
function readtrees2CF(trees::Vector{HybridNetwork};
        quartetfile="none"::AbstractString, whichQ="all"::AbstractString, numQ=0::Integer,
        writeTab=true::Bool, CFfile="none"::AbstractString,
        taxa::AbstractVector=tiplabels(trees),
        writeQ=false::Bool, writeSummary=true::Bool)
    whichQ == "all" || whichQ == "rand" ||
        error("whichQ should be all or rand, not $(whichQ)")
    if(quartetfile == "none")
        readInputData(trees, Symbol(whichQ), numQ, taxa, writeTab, CFfile, writeQ, writeSummary)
    else
        readInputData(trees, quartetfile, Symbol(whichQ), numQ, writeTab, CFfile, writeQ, writeSummary)
    end
end

# ---------------------- descriptive stat for input data ----------------------------------

# function to check how taxa is represented in the input trees
function taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet},s::IO)
    taxaT = tiplabels(trees)
    taxaQ = tiplabels(quartets)
    dif = symdiff(taxaT,taxaQ)
    isempty(dif) ? write(s,"\n same taxa in gene trees and quartets: $(taxaT)\n") :
                   write(s,"\n $(length(dif)) different taxa found in gene trees and quartets. \n Taxa $(intersect(taxaT,dif)) in trees, not in quartets; and taxa $(intersect(taxaQ,dif)) in quartets, not in trees\n")
    u = union(taxaT,taxaQ)
    for taxon in u
        numT = taxonTrees(taxon,trees)
        #numQ = taxonQuartets(taxon,quartets)
        write(s,"Taxon $(taxon) appears in $(numT) input trees ($(round(100*numT/length(trees), digits=2)) %)\n")  #and $(numQ) quartets ($(round(100*numQ/length(quartets), digits=2)) %)\n")
    end
end
taxaTreesQuartets(trees::Vector{HybridNetwork}, quartets::Vector{Quartet}) = taxaTreesQuartets(trees, quartets, stdout)


# function that counts the number of trees in which taxon appears
function taxonTrees(taxon::AbstractString, trees::Vector{HybridNetwork})
    suma = 0
    for t in trees
        suma += in(taxon,t.names) ? 1 : 0
    end
    return suma
end


# function that counts the number of quartets in which taxon appears
function taxonQuartets(taxon::AbstractString, quartets::Vector{Quartet})
    suma = 0
    for q in quartets
        suma += in(taxon,q.taxon) ? 1 : 0
    end
    return suma
end


# function to create descriptive stat from input data, will save in stream sout
# which can be a file or stdout
# default: send to stdout
# pc: only 4-taxon subsets with percentage of gene trees less than pc will be printed (default 70%)
function descData(d::DataCF, sout::IO, pc::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if !isempty(d.tree)
        print(sout,"data consists of $(d.numTrees) gene trees and $(d.numQuartets) 4-taxon subsets\n")
        taxaTreesQuartets(d.tree,d.quartet,sout)
        print(sout,"----------------------------\n\n")
        print(sout,"will print below only the 4-taxon subsets with data from <= $(round((pc)*100, digits=2))% genes\n")
        for q in d.quartet
            percent  = q.ngenes == -1.0 ? 0.0 : round(q.ngenes/d.numTrees*100, digits=2)
            if percent < pc
                print(sout,"4-taxon subset $(q.taxon) obsCF constructed with $(round(q.ngenes)) gene trees ($(percent)%)\n")
            end
        end
        print(sout,"----------------------------\n\n")
    else
        if !isempty(d.quartet)
            print(sout,"data consists of $(d.numQuartets) 4-taxon subsets")
            taxa=tiplabels(d.quartet)
            print(sout,"\nTaxa: $(taxa)\n")
            print(sout,"Number of Taxa: $(length(taxa))\n")
            numQ = binomial(length(taxa),4);
            print(sout,"Maximum number of 4-taxon subsets: $(numQ). Thus, $(round(100*d.numQuartets/numQ, digits=2)) percent of 4-taxon subsets sampled\n")
        end
    end
end

function descData(d::DataCF, filename::AbstractString,pc::Float64)
    println("descriptive stat of input data printed to file $(filename)")
    s = open(filename, "w")
    descData(d,s,pc)
    close(s)
end
descData(d::DataCF, sout::IO=stdout) = descData(d, sout,0.7)
descData(d::DataCF,pc::Float64) = descData(d, stdout,pc)
descData(d::DataCF, filename::AbstractString) = descData(d, filename,0.7)


"""
`summarizedataCF(d::DataCF)`

function to summarize the information contained in a DataCF object. It has the following optional arguments:
- `filename`: if provided, the summary will be saved in the filename, not to screen
- `pc` (number between (0,1)): threshold of percentage of missing genes to identify 4-taxon subsets with fewer genes than the threshold
"""
function summarizedataCF(d::DataCF; filename="none"::AbstractString, pc=0.7::Float64)
    0<=pc<=1 || error("percentage of missing genes should be between 0,1, not: $(pc)")
    if filename == "none"
        descData(d,stdout,pc)
    else
        descData(d,filename,pc)
    end
end


# -------- branch length estimate in coalescent units on species tree ------

"""
    updateBL!(net::HybridNetwork, d::DataCF)

Update internal branch lengths of `net` based on the average quartet concordance
factor (CF) across all quartets that exactly correspond to a given branch:
new branch length = `-log(3/2(1-mean(CF observed in d)))`.
`net` is assumed to be a tree, such that the above equation holds.
"""
function updateBL!(net::HybridNetwork,d::DataCF)
    if !isTree(net)
        @error "updateBL! was created for a tree, and net here is not a tree, so no branch lengths updated"
    end
    parts = edgesParts(net)
    df = makeTable(net,parts,d)
    x = combine(groupby(df, :edge), nrow => :Nquartets,
                :CF => (x -> -log(3/2*(1. - mean(x)))) => :edgeL)
    edges = x[!,:edge]
    lengths = x[!,:edgeL]
    for i in 1:length(edges)
        ind = getIndexEdge(edges[i],net) # helpful error if not found
        if net.edge[ind].length < 0.0 || net.edge[ind].length==1.0
            # readnewicklevel1 changes missing branch length to 1.0
            setLength!(net.edge[ind], (lengths[i] > 0 ? lengths[i] : 0.0))
        end
    end
    for e in net.edge
        if e.length < 0.0 # some edges might have *no* quartet in the data
            setLength!(e, 1.0)
        end
    end
    return x
end


# function to get part1,part2,part3,part4 for each edge in net.edge
# returns a EdgeParts object
function edgesParts(net::HybridNetwork)
    parts = EdgeParts[] #vector to hold part1,...,part4 for each edge
    for e in net.edge
        if(isInternalEdge(e))
            length(e.node) == 2 || error("strange edge with $(length(e.node)) nodes instead of 2")
            n1 = e.node[1]
            n2 = e.node[2]
            e11,e12 = hybridEdges(n1,e)
            e21,e22 = hybridEdges(n2,e)
            part1 = Node[]
            part2 = Node[]
            part3 = Node[]
            part4 = Node[]
            getDescendants!(getOtherNode(e11,n1),e11,part1)
            getDescendants!(getOtherNode(e12,n1),e12,part2)
            getDescendants!(getOtherNode(e21,n2),e21,part3)
            getDescendants!(getOtherNode(e22,n2),e22,part4)
            push!(parts, EdgeParts(e.number,part1,part2,part3,part4))
        end
    end
    return parts
end

# Traverse the network from a node towards an edge, following major edges only.
# The other `getDescendants` updates the nodes `.inCycle`, to update partitions later.
function getDescendants!(
    node::Node,
    edge::Edge,
    descendants::Array{Node,1}
)
    if node.leaf
        push!(descendants, node)
    else
        for e in node.edge
            if !isEqual(edge,e) && e.ismajor
                other = getOtherNode(e,node);
                getDescendants!(other,e, descendants);
            end
        end
    end
end


# similar to `getDescendants!` above, but uses a vector of edges instead of nodes
# finds the partition corresponding to the node and edge in the cycle
# used in chooseEdgesGamma and to set net.partition
# cycleNum is a variable that will save another hybrid node number if found
function getDescendants!(
    node::Node,
    edge::Edge,
    descendants::Vector{Edge},
    cycleNum::Vector{Int}
)
    @debug "getDescendants of node $(node.number) and edge $(edge.number)"
    if inCycle(node) != -1
        push!(cycleNum,inCycle(node))
    elseif !node.leaf && inCycle(node) == -1
        for e in node.edge
            if(!isEqual(edge,e) && e.ismajor)
                push!(descendants,e)
                getDescendants!(getOtherNode(e,node),e,descendants,cycleNum)
            end
        end
    end
end


# function to make table to later use in updateBL
# uses vector parts obtained from edgeParts function
function makeTable(net::HybridNetwork, parts::Vector{EdgeParts},d::DataCF)
    df = DataFrame(edge=Int[],t1=AbstractString[],t2=AbstractString[],t3=AbstractString[],t4=AbstractString[],resolution=AbstractString[],CF=Float64[])
    sortedDataQ = [sort(q.taxon) for q in d.quartet]
    for p in parts #go over internal edges too
        for t1 in p.part1
            for t2 in p.part2
                for t3 in p.part3
                    for t4 in p.part4
                        tx1 = net.names[t1.number]
                        tx2 = net.names[t2.number]
                        tx3 = net.names[t3.number]
                        tx4 = net.names[t4.number]
                        nam = [tx1,tx2,tx3,tx4]
                        snam = sort(nam)
                        row = findall(isequal(snam), sortedDataQ)
                        for r in row # nothing if tax set not found: length(row)=0
                          col,res = resolution(nam,d.quartet[r].taxon)
                          push!(df, [p.edgenum,tx1,tx2,tx3,tx4,res,d.quartet[r].obsCF[col]])
                        end
                    end
                end
            end
        end
    end
    return df
end


# function to determine the resolution of taxa picked from part1,2,3,4 and DataCF
# names: taxa from part1,2,3,4
# rownames: taxa from table of obsCF
function resolution(names::Vector{<:AbstractString},rownames::Vector{<:AbstractString})
    length(names) == length(rownames) || error("names and rownames should have the same length")
    length(names) == 4 || error("names should have 4 entries, not $(length(names))")
    bin = [n == names[1] || n == names[2] ? 1 : 0 for n in rownames]
    if(bin == [1,1,0,0] || bin == [0,0,1,1])
        return 1,"12|34"
    elseif(bin == [1,0,1,0] || bin == [0,1,0,1])
        return 2,"13|24"
    elseif(bin == [1,0,0,1] || bin == [0,1,1,0])
        return 3,"14|23"
    else
        error("strange resolution $(bin)")
    end
end


# function to extract a quartet from a matrix M
# obtained from tree2Matrix (defined in file compareNetworks.jl)
# this function is meant to replace extractQuartet! in calculateObsCFAll
# input: Quartet, Matrix, vector of taxa names
# returns 1 if quartet found is 12|34, 2 if 13|24, 3 if 14|23, and 0 if not found
function extractQuartetTree(q::Quartet, M::Matrix{Int},S::Union{Vector{<:AbstractString},Vector{Int}})
    @debug "extractQuartet: $(q.taxon)"
    @debug "matrix: $(M)"
    inds = indexin(q.taxon, S)
    if any(isnothing, inds)
        error("some taxon in quartet $(q.taxon) not found in list of all species $(S)")
    end
    subM = M[:, inds.+1]
    @debug "subM: $(subM)"
    for r in 1:size(subM,1) #rows in subM
        @debug "subM[r,:]: $(subM[r,:])"
        if subM[r,:] == [0,0,1,1] || subM[r,:] == [1,1,0,0]
            return 1
        elseif subM[r,:] == [0,1,0,1] || subM[r,:] == [1,0,1,0]
            return 2
        elseif subM[r,:] == [0,1,1,0] || subM[r,:] == [1,0,0,1]
            return 3
        end
    end
    return 0
end


# function that will give the qth quartet without making a list of all quartets
# input: n number of taxa, q desired index of quartet
# returns vector of int, e.g. 1234
function whichQuartet(n::Int, q::Int)
    p = 4
    q <= binom(n,p) || error("the index for the quartet $(q) needs to be less than choose(n,4)=$(binom(n,p))")
    n > 4 || error("there must be at least 5 taxa, not $(n)")
    quartet = Int[]
    while(n > 1)
        abs = binom(n-1,p) #fixit: we don't want to compute this, we want to look for it in a table
        if(q > abs)
            push!(quartet,n)
            n -= 1
            p -= 1
            q = q-abs
        else
            n -= 1
        end
    end
    if(length(quartet) == 3)
        push!(quartet,1)
    end
    quartet = quartet[[4,3,2,1]] #sort
    return quartet
end


# function to write a quartet on integer to taxon names
# it creates a Quartet type
# input: list of taxa names, vector of integers (each integer corresponds to a taxon name),
# num is the number of the Quartet
# assumes taxa is sorted already
function createQuartet(taxa::Union{Vector{<:AbstractString},Vector{Int}}, qvec::Vector{Int}, num::Int)
    length(qvec) == 4 || error("a quartet should have only 4 taxa, not $(length(qvec))")
    names = String[]
    for i in qvec
        i <= length(taxa) || error("want taxon number $(i) in list of taxon names $(taxa) which has only $(length(taxa)) names")
        push!(names,string(taxa[i]))
    end
    return Quartet(num,names,[1.0,0.0,0.0])
end
