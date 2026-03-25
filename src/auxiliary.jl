# auxiliary functions for all the other methods
# originally in functions.jl
# Claudia February 2015
#####################
# ----- aux general functions ---------------

function isInternalEdge(edge::Edge)
    length(edge.node) == 2 || error("edge $(edge.number) has $(length(edge.node)) nodes, should be 2")
    return !edge.node[1].leaf && !edge.node[2].leaf
end


# ----------------------------------------------------------------------------------------

# setLength
# warning: allows to change edge length for istIdentifiable=false
#          but issues a warning
# negative=true means it allows negative branch lengths
function setLength!(edge::Edge, new_length::Number, negative::Bool)
    (negative || new_length >= 0) || error("length has to be nonnegative: $(new_length), cannot set to edge $(edge.number)")
    new_length >= -0.4054651081081644 || error("length can be negative, but not too negative (greater than -log(1.5)) or majorCF<0: new length is $(new_length)")
    #println("setting length $(new_length) to edge $(edge.number)")
    if(new_length > 10.0)
        new_length = 10.0;
    end
    edge.length = new_length;
    edge.y = exp(-new_length);
    edge.z = 1.0 - edge.y;
    #istIdentifiable(edge) || @warn "set edge length for edge $(edge.number) that is not identifiable"
    return nothing
end

"""
    setLength!(edge, newlength)`

Set the length of `edge`, and set `edge.y` and `edge.z` accordingly.
Warning: specific to `SNaQ.jl`.
Consider [`PhyloNetworks.setlengths!`](@extref) from `PhyloNetworks` for a more generic tool.

- The new length is censored to 10: if the new length is above 10,
  the edge's length will be set to 10. Lengths are interpreted in coalescent
  units, and 10 is close to infinity: near perfect gene tree concordance.
  10 is used as an upper limit to coalescent units that can be reliably estimated.
- The new length is allowed to be negative, but must be greater than -log(1.5),
  to ensure that the major quartet concordance factor (1 - 2/3 exp(-length)) is >= 0.
"""
setLength!(edge::Edge, new_length::Number) = setLength!(edge, new_length, false)



"""
    setGammaBLfromGammaz!(node, network)

Update the γ values of the two sister hybrid edges in a bad diamond I, given the `gammaz` values
of their parent nodes, and update the branch lengths t1 and t2 of their parent edges
(those across from the hybrid nodes), in such a way that t1=t2 and that these branch lengths
and γ values are consistent with the `gammaz` values in the network.

Similar to the first section of [`undoGammaz!`](@ref),
but does not update anything else than γ and t's.
Unlike `undoGammaz!`, no error if non-hybrid `node` or not at bad diamond I.
"""
function setGammaBLfromGammaz!(node::Node, net::HybridNetwork)
    if !isBadDiamondI(node) || !node.hybrid
        return nothing
    end
    edge_maj, edge_min, tree_edge2 = hybridEdges(node);
    other_maj = getOtherNode(edge_maj,node);
    other_min = getOtherNode(edge_min,node);
    edgebla,tree_edge_incycle1,tree_edge = hybridEdges(other_min);
    edgebla,tree_edge_incycle2,tree_edge = hybridEdges(other_maj);
    if(approxEq(gammaz(other_maj),0.0) && approxEq(gammaz(other_min),0.0))
        edge_maj.gamma = 1.0 # γ and t could be anything if both gammaz are 0
        edge_min.gamma = 0.0 # will set t's to 0 and minor γ to 0.
        newt = 0.0
    else
        ((approxEq(gammaz(other_min),0.0) || gammaz(other_min) >= 0.0) &&
         (approxEq(gammaz(other_maj),0.0) || gammaz(other_maj) >= 0.0)    ) ||
            error("bad diamond I in node $(node.number) but missing (or <0) gammaz")
        ztotal = gammaz(other_maj) + gammaz(other_min)
        edge_maj.gamma = gammaz(other_maj) / ztotal
        edge_min.gamma = gammaz(other_min) / ztotal
        newt = -log(1-ztotal)
    end
    setLength!(tree_edge_incycle1,newt)
    setLength!(tree_edge_incycle2,newt)
end

# function to find if a given partition is in net.partition
function isPartitionInNet(net::HybridNetwork,desc::Vector{Edge},cycle::Vector{Int})
    for p in net.partition
        if(sort(cycle) == sort(p.cycle))
            if(sort([e.number for e in desc]) == sort([e.number for e in p.edges]))
                return true
            end
        end
    end
    return false
end

# function to check if a partition is already in net.partition
# used in updatePartition
function isPartitionInNet(net::HybridNetwork,partition::Partition)
    if(isempty(net.partition))
        return false
    end
    for p in net.partition
        cycle = isempty(setdiff(p.cycle,partition.cycle)) && isempty(setdiff(partition.cycle,p.cycle))
        edges = isempty(setdiff([n.number for n in p.edges],[n.number for n in partition.edges])) && isempty(setdiff([n.number for n in partition.edges],[n.number for n in p.edges]))
        if(cycle && edges)
            return true
        end
    end
    return false
end

"""
    sorttaxa!(DataFrame, columns)

Reorder the 4 taxa and reorders the observed concordance factors accordingly, on each row of
the data frame. If `columns` is ommitted, taxon names are assumed to be in columns 1-4 and
CFs are assumed to be in columns 5-6 with quartets in this order: `12_34`, `13_24`, `14_23`.
Does **not** reorder credibility interval values, if present.

    sorttaxa!(DataCF)
    sorttaxa!(Quartet, permutation_tax, permutation_cf)

Reorder the 4 taxa in each element of the DataCF `quartet`. For a given Quartet,
reorder the 4 taxa in its fields `taxon` (if non-empty)
and reorder the 3 concordance values accordingly, in `obsCF`

`permutation_tax` and `permutation_cf` should be vectors of short integers (Int8) of length 4 and 3
respectively, whose memory allocation gets reused. Their length is *not checked*.
"""
function sorttaxa!(dat::DataCF)
    ptax = Array{Int8}(undef, 4) # to hold the sort permutations
    pCF  = Array{Int8}(undef, 3)
    for q in dat.quartet
        sorttaxa!(q, ptax, pCF)
    end
end

function sorttaxa!(df::DataFrame, co=Int[]::Vector{Int})
    if length(co)==0
        co = collect(1:7)
    end
    length(co) > 6 || error("column vector must be of length 7 or more")
    ptax = Array{Int8}(undef, 4)
    pCF  = Array{Int8}(undef, 3)
    taxnam = Array{eltype(df[!,co[1]])}(undef, 4)
    for i in axes(df,1)
        for j=1:4 taxnam[j] = df[i,co[j]]; end
        sortperm!(ptax, taxnam)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF according to taxon permutation
        df[i,co[1]], df[i,co[2]], df[i,co[3]], df[i,co[4]] = taxnam[ptax[1]], taxnam[ptax[2]], taxnam[ptax[3]], taxnam[ptax[4]]
        df[i,co[5]], df[i,co[6]], df[i,co[7]] = df[i,co[pCF[1]+4]], df[i,co[pCF[2]+4]], df[i,co[pCF[3]+4]]
    end
    return df
end

function sorttaxa!(qua::Quartet, ptax::Vector{Int8}, pCF::Vector{Int8})
    qt = qua.taxon
    if length(qt)==4
        sortperm!(ptax, qt)
        sorttaxaCFperm!(pCF, ptax) # update permutation pCF accordingly
        qt[1], qt[2], qt[3], qt[4] = qt[ptax[1]], qt[ptax[2]], qt[ptax[3]], qt[ptax[4]]
        qua.obsCF[1], qua.obsCF[2], qua.obsCF[3] = qua.obsCF[pCF[1]], qua.obsCF[pCF[2]], qua.obsCF[pCF[3]]
    elseif length(qt)!=0
        error("Quartet with $(length(qt)) taxa")
    end
    return qua
end

# find permutation pCF of the 3 CF values: 12_34, 13_24, 14_23. 3!=6 possible permutations
# ptax = one of 4!=24 possible permutations on the 4 taxon names
# kernel: pCF = identity if ptax = 1234, 2143, 3412 or 4321
# very long code, but to minimize equality checks at run time
function sorttaxaCFperm!(pcf::Vector{Int8}, ptax::Vector{Int8})
    if ptax[1]==1
        if     ptax[2]==2
            pcf[1]=1
            if  ptax[3]==3 # ptax = 1,2,3,4
                pcf[2]=2; pcf[3]=3
            else           # ptax = 1,2,4,3
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==3
            pcf[1]=2
            if  ptax[3]==2 # ptax = 1,3,2,4
                pcf[2]=1; pcf[3]=3
            else           # ptax = 1,3,4,2
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==4
            pcf[1]=3
            if  ptax[3]==2 # ptax = 1,4,2,3
                pcf[2]=1; pcf[3]=2
            else           # ptax = 1,4,3,2
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==2
        if     ptax[2]==1
            pcf[1]=1
            if  ptax[3]==4 # ptax = 2,1,4,3
                pcf[2]=2; pcf[3]=3
            else           # ptax = 2,1,3,4
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==4
            pcf[1]=2
            if  ptax[3]==1 # ptax = 2,4,1,3
                pcf[2]=1; pcf[3]=3
            else           # ptax = 2,4,3,1
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==3
            pcf[1]=3
            if  ptax[3]==1 # ptax = 2,3,1,4
                pcf[2]=1; pcf[3]=2
            else           # ptax = 2,3,4,1
                pcf[2]=2; pcf[3]=1
            end
        end
    elseif ptax[1]==3
        if     ptax[2]==4
            pcf[1]=1
            if  ptax[3]==1 # ptax = 3,4,1,2
                pcf[2]=2; pcf[3]=3
            else           # ptax = 3,4,2,1
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==1
            pcf[1]=2
            if  ptax[3]==4 # ptax = 3,1,4,2
                pcf[2]=1; pcf[3]=3
            else           # ptax = 3,1,2,4
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==2
            pcf[1]=3
            if  ptax[3]==4 # ptax = 3,2,4,1
                pcf[2]=1; pcf[3]=2
            else           # ptax = 3,2,1,4
                pcf[2]=2; pcf[3]=1
            end
        end
    else # ptax[1]==4
        if     ptax[2]==3
            pcf[1]=1
            if  ptax[3]==2 # ptax = 4,3,2,1
                pcf[2]=2; pcf[3]=3
            else           # ptax = 4,3,1,2
                pcf[2]=3; pcf[3]=2
            end
        elseif ptax[2]==2
            pcf[1]=2
            if  ptax[3]==3 # ptax = 4,2,3,1
                pcf[2]=1; pcf[3]=3
            else           # ptax = 4,2,1,3
                pcf[2]=3; pcf[3]=1
            end
        else # ptax[2]==1
            pcf[1]=3
            if  ptax[3]==3 # ptax = 4,1,3,2
                pcf[2]=1; pcf[3]=2
            else           # ptax = 4,1,2,3
                pcf[2]=2; pcf[3]=1
            end
        end
    end
end

# function to check in an edge is in an array by comparing
# edge numbers (could use isEqual for adding comparisons of gammaz and inCycle)
# needed for updateHasEdge
function isEdgeNumIn(edge::Edge,array::Array{Edge,1})
    enum = edge.number
    return any(e -> e.number == enum, array)
end

# function to check in a leaf is in an array by comparing
# the numbers (uses isEqual)
# needed for updateHasEdge
function isNodeNumIn(node::Node,array::Array{Node,1})
    return all((e->!isEqual(node,e)), array) ? false : true
end



#------------------------------------
function citation()
    bibfile = joinpath(@__DIR__, "..", "CITATION.bib")
    out = readlines(bibfile)
    println("Bibliography in bibtex format also in CITATION.bib")
    println(join(out,'\n'))
end
