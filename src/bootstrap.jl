# julia functions for bootstrap
# Claudia October 2015
# Cecile April 2016



"""
    sampleCFfromCI(data frame, seed=0)
    sampleCFfromCI!(data frame, seed=0)

Read a data frame containing CFs and their credibility intervals, and
sample new obsCF uniformly within the CIs.
These CFs are then rescaled to sum up to 1 for each 4-taxon sets.
Return a data frame with taxon labels in first 4 columns, sampled obs CFs in columns 5-7
and credibility intervals in columns 8-13.

- The non-modifying function creates a new data frame (with re-ordered columns) and returns it.
  If `seed=-1`, the new df is a deep copy of the input df, with no call to the random
  number generator. Otherwise, `seed` is passed to the modifying function.
- The modifying function overwrites the input data frame with the sampled CFs and returns it.
  If `seed=0`, the random generator is seeded from the clock. Otherwise the random generator
  is seeded using `seed`.

Warning: the modifying version does *not* check the data frame: assumes correct columns.

optional argument: `delim=','` by default: how columns are delimited.
"""
function sampleCFfromCI(df::DataFrame, seed=0::Integer)
    @debug "order of columns should be: t1,t2,t3,t4,cf1234,cf1324,cf1423,cf1234LO,cf1234HI,..."
    size(df,2) == 13 || size(df,2) == 14 || @warn "sampleCFfromCI function assumes table from TICR: CF, CFlo, CFhi"
    obsCFcol = [findfirst(isequal(:CF12_34), DataFrames.propertynames(df)),
                findfirst(isequal(:CF13_24), DataFrames.propertynames(df)),
                findfirst(isequal(:CF14_23), DataFrames.propertynames(df))]
    nothing ∉ obsCFcol || error("""CF columns were not found: should be named like 'CF12_34'""")
    obsCFcol == [5,8,11] ||
        @warn """CF columns were found, but not in the expected columns.
                Lower/upper bounds of credibility intervals assumed in columns 6,7, 9,10 and 12,13."""
    colsTa = [1,2,3,4]        # column numbers for taxon names
    colsCI = [6,7,9,10,12,13] # for lower/upper CI bounds
    length(findall(in(obsCFcol), colsTa)) ==0 ||
        error("CFs found in columns 1-4 where taxon labels are expected")
    length(findall(in(obsCFcol), colsCI)) ==0 ||
        error("CFs found in columns where credibility intervals are expected")
    newdf = df[:, [colsTa; obsCFcol; colsCI] ]
    if seed==-1
      return newdf
    else
      return sampleCFfromCI!(newdf::DataFrame, seed)
    end
end

function sampleCFfromCI!(df::DataFrame, seed=0::Integer)
    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
        println("using seed $(seed) for bootstrap table")
    end
    Random.seed!(seed)
    for i in axes(df,1)
        c1 = (df[i, 9]-df[i, 8])*rand()+df[i, 8]
        c2 = (df[i,11]-df[i,10])*rand()+df[i,10]
        c3 = (df[i,13]-df[i,12])*rand()+df[i,12]
        suma = c1+c2+c3
        df[i,5] = c1/suma
        df[i,6] = c2/suma
        df[i,7] = c3/suma
    end
    return df
end

sampleCFfromCI(file::AbstractString; delim=','::Char,seed=0::Integer) =
    sampleCFfromCI(DataFrame(CSV.File(file, delim=delim); copycols=false),seed)

# function that will do bootstrap of snaq estimation in series
# it repeats optTopRuns nrep times
# it has the same arguments as optTopRuns except for:
# - need data table of CF with conf intervals (instead of d DataCF),
#   or vector of vector of HybridNetworks
# - new argument nrep: number of bootstrap replicates (default 10)
# - new argument runs2: percentage of bootstrap replicates to start in the best network, by default 0.25
# - new argument bestNet: to start the optimization. if prcnet>0.0 and bestNet is not input as argument from a previous run, it will estimate it inside
# - quartetfile if it was used in original data ("none" if all quartets used)
# recall: optTopRuns! does *not* modify its input starting network
function optTopRunsBoot(currT0::HybridNetwork, data::Union{DataFrame,Vector{Vector{HybridNetwork}}},
                        hmax::Integer, liktolAbs::Float64, Nfail::Integer, ftolRel::Float64,ftolAbs::Float64,xtolRel::Float64,xtolAbs::Float64,
                        verbose::Bool, closeN::Bool, Nmov0::Vector{Int},
                        runs1::Integer, outgroup::AbstractString, filename::AbstractString, seed::Integer, probST::Float64,
                        nrep::Integer, runs2::Integer, bestNet::HybridNetwork, quartetfile::AbstractString,
                        probQR::AbstractFloat, propQuartets::AbstractFloat)
    println("BOOTSTRAP OF SNAQ ESTIMATION")
    writelog = true
    if filename != ""
        logfile = open(string(filename,".log"),"w")
        write(logfile, "BOOTSTRAP OF SNAQ ESTIMATION \n")
    else
        writelog = false
        logfile = stdout
    end

    inputastrees = isa(data, Vector{Vector{HybridNetwork}})
    inputastrees || isa(data, DataFrame) ||
        error("Input data not recognized: $(typeof(data))")

    if runs1>0 && runs2>0
        str = """Will use this network as starting topology for $runs1 run(s) for each bootstrap replicate:
                 $(writenewick_level1(currT0))
                 and this other network for $runs2 run(s):
                 $(writenewick_level1(bestNet))
                 """
        writelog && write(logfile, str)
        print(str)
    end

    if inputastrees # allocate memory, to be re-used later
        newtrees = samplebootstrap_multiloci(data, row=1)
        newd = readtrees2CF(newtrees, quartetfile=quartetfile, writeTab=false, writeSummary=false)
        taxa = tiplabels(newtrees)
    else
        newdf = sampleCFfromCI(data, -1) # column names check, newdf has obsCF in columns 5-7
        # seed=-1: deep copy only, no rand()
        newd = readtableCF!(newdf, collect(1:7)) # allocate memory for DataCF object
    end
    if runs1>0 && isTree(currT0) # get rough first estimate of branch lengths in startnet
        updateBL!(currT0, newd)
    end

    if seed == 0
        t = time()/1e9
        a = split(string(t),".")
        seed = parse(Int,a[2][end-4:end]) #better seed based on clock
    end
    println("main seed $(seed)")
    writelog && write(logfile,"\nmain seed $(seed)\n")
    writelog && flush(logfile)
    Random.seed!(seed)
    seedsData = round.(Int,floor.(rand(nrep)*100000)) # seeds to sample bootstrap data
    if runs1>0
        seeds = round.(Int,floor.(rand(nrep)*100000)) # seeds for all optimizations from currT0
    end
    if runs2>0
      seedsOtherNet = round.(Int,floor.(rand(nrep)*100000)) # for runs starting from other net
    end

    bootNet = HybridNetwork[]

    writelog && write(logfile,"\nBEGIN: $(nrep) replicates\n$(Libc.strftime(time()))\n")
    writelog && flush(logfile)
    for i in 1:nrep
        str = "\nbegin replicate $(i)\nbootstrap data simulation: seed $(seedsData[i])\n"
        writelog && write(logfile, str)
        print(str)
        if !inputastrees
            sampleCFfromCI!(newdf, seedsData[i])
            readtableCF!(newd, newdf, [5,6,7])
        else
            samplebootstrap_multiloci!(newtrees, data, seed=seedsData[i])
            calculateObsCFAll!(newd,taxa) # do not use readtrees2CF: to save memory and gc time
        end
        if runs1>0
            str = "estimation, $runs1 run" * (runs1>1 ? "s" : "") * ": seed $(seeds[i])\n"
            writelog && write(logfile, str)
            print(str)
            rootname = ""
            @debug begin rootname = string(filename,"_",i);
                         "rootname set to $rootname"; end
            net1 = optTopRuns!(currT0, liktolAbs, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs1, outgroup,
                               rootname,seeds[i],probST,probQR,propQuartets)
            if runs2==0
                net = net1
            end
        end
        if runs2>0
            str = "estimation, $runs2 run" * (runs2>1 ? "s" : "") * " starting from other net: seed $(seedsOtherNet[i])\n"
            writelog && write(logfile, str)
            print(str)
            rootname = ""
            @debug begin rootname = string(filename,"_",i,"_startNet2");
                         "rootname set to $rootname"; end
            net2 = optTopRuns!(bestNet, liktolAbs, Nfail, newd, hmax,ftolRel, ftolAbs, xtolRel, xtolAbs, verbose, closeN, Nmov0, runs2, outgroup,
                               rootname,seedsOtherNet[i],probST,probQR,propQuartets)
            if runs1==0
                net = net2
            end
        end
        if runs1>0 && runs2>0
            net = (loglik(net1) < loglik(net2) ? net1 : net2)
        end
        writelog && flush(logfile)

        push!(bootNet, net)
        str = (outgroup=="none" ? writenewick_level1(net) : writenewick_level1(net,outgroup))
        if writelog
            write(logfile, str)
            write(logfile,"\n")
            flush(logfile)
        end
        println(str) # net also printed by each optTopRuns! but separately from the 2 starting points
    end # of the nrep bootstrap replicates
    writelog && close(logfile)

    if writelog
      s = open(string(filename,".out"),"w")
      for n in bootNet
        if outgroup == "none"
            write(s,"$(writenewick_level1(n))\n")
        else
            write(s,"$(writenewick_level1(n,outgroup))\n")
        end
        # "with -loglik $(loglik(n))" not printed: not comparable across bootstrap networks
      end
      close(s)
    end
    return bootNet
end

# like snaq, only calls optTopRunsBoot
# undocumented arguments: closeN, Nmov0
"""
    bootsnaq(T::HybridNetwork, df::DataFrame)
    bootsnaq(T::HybridNetwork, vector of tree lists)

Bootstrap analysis for SNaQ.
Bootstrap data can be quartet concordance factors (CF),
drawn from sampling uniformly in their credibility intervals,
as given in the data frame `df`.
Alternatively, bootstrap data can be gene trees sampled from
a vector of tree lists: one list of bootstrap trees per locus
(see [`readmultinewick_files`](https://juliaphylo.github.io/PhyloNetworks.jl/stable/lib/public/#PhyloNetworks.readmultinewick_files-Tuple{AbstractString}) in `PhyloNetworks` to generate this,
from a file containing a list of bootstrap files: one per locus).

From each bootstrap replicate, a network
is estimated with [`snaq!`](@ref), with a search starting from topology `T`.
Optional arguments include the following, with default values in parentheses:

- `hmax` (1): max number of reticulations in the estimated networks
- `nrep` (10): number of bootstrap replicates.
- `runs` (10): number of independent optimization runs for each replicate
- `filename` ("bootsnaq"): root name for output files. No output files if "".
- `seed` (0 to get a random seed from the clock): seed for random number generator
- `otherNet` (empty): another starting topology so that each replicate will start prcnet% runs on otherNet and (1-prcnet)% runs on `T`
- `prcnet` (0): percentage of runs starting on `otherNet`; error if different than 0.0, and otherNet not specified.
- `ftolRel`, `ftolAbs`, `xtolRel`, `xtolAbs`, `liktolAbs`, `Nfail`,
  `probST`, `verbose`, `outgroup`: see `snaq!`, same defaults.

If `T` is a tree, its branch lengths are first optimized roughly with [`updateBL!`](@ref)
(by using the average CF of all quartets defining each branch and calculating the coalescent units
corresponding to this quartet CF).
If `T` has one or more reticulations, its branch lengths are taken as is to start the search.
The branch lengths of `otherNet` are always taken as is to start the search.
"""
function bootsnaq(startnet::HybridNetwork, data::Union{DataFrame,Vector{Vector{HybridNetwork}}};
                  hmax=1::Integer, liktolAbs=likAbs::Float64, Nfail=numFails::Integer,
                  ftolRel=fRel::Float64, ftolAbs=fAbs::Float64, xtolRel=xRel::Float64, xtolAbs=xAbs::Float64,
                  verbose=false::Bool, closeN=true::Bool, Nmov0=numMoves::Vector{Int},
                  runs=10::Integer, outgroup="none"::AbstractString, filename="bootsnaq"::AbstractString,
                  seed=0::Integer, probST=0.3::Float64, nrep=10::Integer, prcnet=0.0::Float64,
                  otherNet=HybridNetwork()::HybridNetwork, quartetfile="none"::AbstractString,
                  probQR::AbstractFloat=0.0, propQuartets::AbstractFloat=1.0)

    inputastrees = isa(data, Vector{Vector{HybridNetwork}})
    inputastrees || isa(data, DataFrame) ||
        error("Input data not recognized: $(typeof(data))")

    if !inputastrees
    (DataFrames.propertynames(data)[[6,7,9,10,12,13]] == [:CF12_34_lo,:CF12_34_hi,:CF13_24_lo,:CF13_24_hi,:CF14_23_lo,:CF14_23_hi]) ||
      @warn """assume table with CI from TICR: CFlo, CFhi in columns 6,7; 9,10; and 12,13.
              Found different column names: $(DataFrames.names(data)[[6,7,9,10,12,13]])"""
    else # check 1+ genes, each with 1+ trees, all with h=0.
        ngenes = length(data)
        ngenes > 0 || error("empty list of bootstrap trees (0 genes)")
        for igene in 1:ngenes
            btr = data[igene]
            length(btr) > 0 || error("no bootstrap trees for $(igene)th gene")
            for itree in eachindex(btr)
                btr[itree].numhybrids == 0 || error("network $itree is not a tree for $(igene)th gene")
            end
        end
    end
    prcnet >= 0 || error("percentage of times to use the best network as starting topology should be positive: $(prcnet)")
    prcnet = (prcnet <= 1.0) ? prcnet : prcnet/100
    runs2 = round(Int, runs*prcnet)            # runs starting from otherNet
    runs1 = runs - runs2                       # runs starting from startnet

    if runs1>0
        startnet=readnewicklevel1(writenewick_level1(startnet)) # does not modify startnet outside
        flag = checkNet(startnet,true) # light checking only
        flag && error("starting topology suspected not level-1")
        try
            checkNet(startnet)
        catch err
            println("starting topology is not of level 1:")
            rethrow(err)
        end
    end
    runs2 == 0 || otherNet.numtaxa > 0 ||
        error("""otherNet not given and prcnet>0. Please set prcnet to 0 to start optimizations
                from the same network always, or else provide an other network "otherNet"
                to start the optimization from this other network in pcrnet % of runs.""")
    if runs2 > 0
        otherNet=readnewicklevel1(writenewick_level1(otherNet))
        flag = checkNet(otherNet,true) # light checking only
        flag && error("starting topology 'otherNet' suspected not level-1")
        try
            checkNet(otherNet)
        catch err
            println("starting topology 'otherNet' not a level 1 network:")
            rethrow(err)
        end
    end

    # for multiple alleles: expand into two leaves quartets like sp1 sp1 sp2 sp3.
    if (@isdefined originald) && !isempty(originald.repSpecies) ## not defined if treefile empty, but not needed
        expandLeaves!(originald.repSpecies,startnet)
    end

    optTopRunsBoot(startnet,data,hmax, liktolAbs, Nfail,ftolRel, ftolAbs, xtolRel, xtolAbs,
                   verbose, closeN, Nmov0, runs1, outgroup, filename,
                   seed, probST, nrep, runs2, otherNet, quartetfile, probQR, propQuartets)
end

