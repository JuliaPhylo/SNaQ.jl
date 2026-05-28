import Base: getproperty, getfield, getindex, setproperty!, show


# Edge getters/setters
istIdentifiable(e::Edge) = e.boole1 # true if the parameter t (length) for this edge is identifiable as part of a network
fromBadDiamondI(e::Edge) = e.boole2 # true if the edge came from deleting a bad diamond I hybridization
inCycle(e::Edge) = e.inte1 # = Hybrid node number if this edge is part of a cycle created by such hybrid node; -1 if not part of cycle. used to add new hybrid edge. updated after edge is part of a network

istIdentifiable!(e::Edge, b::Bool) = (e.boole1 = b)
fromBadDiamondI!(e::Edge, b::Bool) = (e.boole2 = b)
inCycle!(e::Edge, i::Int) = (e.inte1 = i)


# Node getters/setters
hasHybEdge(n::Node) = n.booln1 # is there a hybrid edge in edge? only needed when hybrid=false (tree node)
isBadDiamondI(n::Node) = n.booln2 # for hybrid node, is it bad diamond case I, update in updateGammaz!
isBadDiamondII(n::Node) = n.booln3 # for hybrid node, is it bad diamond case II, update in updateGammaz!
isExtBadTriangle(n::Node) = n.booln4 # for hybrid node, is it extremely bad triangle, udpate in updateGammaz!
isVeryBadTriangle(n::Node) = n.booln5 # for hybrid node, is it very bad triangle, udpate in updateGammaz!
isBadTriangle(n::Node) = n.booln6 # for hybrid node, is it bad triangle, udpate in updateGammaz!
k(n::Node) = n.intn2 # num nodes in cycle, only stored in hybrid node, updated after node becomes part of network; default -1
typeHyb(n::Node) = n.int8n3 # type of hybridization (1,2,3,4, or 5), needed for quartet network only. default -1
inCycle(n::Node) = n.intn1 # = hybrid node if this node is part of a cycle created by such hybrid node, -1 if not part of cycle
gammaz(n::Node) = n.fvalue # notes file for explanation. gammaz if tree node, gamma2z if hybrid node; updated after node is part of network with updateGammaz!

hasHybEdge!(n::Node, b::Bool) = (n.booln1 = b)
isBadDiamondI!(n::Node, b::Bool) = (n.booln2 = b)
isBadDiamondII!(n::Node, b::Bool) = (n.booln3 = b)
isExtBadTriangle!(n::Node, b::Bool) = (n.booln4 = b)
isVeryBadTriangle!(n::Node, b::Bool) = (n.booln5 = b)
isBadTriangle!(n::Node, b::Bool) = (n.booln6 = b)
k!(n::Node, i::Int) = (n.intn2 = i)
typeHyb!(n::Node, i::Int) = typeHyb!(n, Int8(i))
typeHyb!(n::Node, i::Int8) = (n.int8n3 = i)
inCycle!(n::Node, i::Int) = (n.intn1 = i)
gammaz!(n::Node, f::Real) = (n.fvalue = f)


# HybridNetwork getters/setters
visited(h::HybridNetwork) = h.vec_bool
edges_changed(h::HybridNetwork) = h.vec_edge
nodes_changed(h::HybridNetwork) = h.vec_node
ht(h::HybridNetwork) = h.vec_float # vector of parameters to optimize
numht(h::HybridNetwork) = h.vec_int2 # vector of number of the hybrid nodes and edges in ht e.g. [3,6,8,...], 2 hybrid nodes 3,6, and edge 8 is the 1st identifiable
numBad(h::HybridNetwork) = h.intg1 # number of bad diamond I hybrid nodes, set as 0
hasVeryBadTriangle(h::HybridNetwork) = h.boolg1 # true if the network has extremely/very bad triangles that should be ignored
index(h::HybridNetwork) = h.vec_int3 #index in net.edge, net.node of elements in net.ht to make updating easy
loglik(h::HybridNetwork) = h.fscore # composite log-likelihood (higher = better fit)
blacklist(h::HybridNetwork) = h.vec_int4 # reusable array of integers, used in afterOptBL
cleaned(h::HybridNetwork) = h.boolg2 # attribute to know if the network has been cleaned after readm default false

visited!(h::HybridNetwork, v::BitVector) = visited!(h, Vector{Bool}(v))
visited!(h::HybridNetwork, v::Array{Bool, 1}) = (h.vec_bool = v)
edges_changed!(h::HybridNetwork, v::Array{Edge, 1}) = (h.vec_edge = v)
nodes_changed!(h::HybridNetwork, v::Array{Node, 1}) = (h.vec_node = v)
ht!(h::HybridNetwork, v::Vector{<:Real}) = (h.vec_float = v)
numht!(h::HybridNetwork, v::Vector{Int}) = (h.vec_int2 = v)
numBad!(h::HybridNetwork, i::Int) = (h.intg1 = i)
hasVeryBadTriangle!(h::HybridNetwork, b::Bool) = (h.boolg1 = b)
index!(h::HybridNetwork, v::Vector{Int}) = (h.vec_int3 = v)

"""
Sets the composite log-likelihood of network `h` to `f`.
Higher values indicate better fit.
"""
loglik!(h::HybridNetwork, f::Real) = (h.fscore = f)
blacklist!(h::HybridNetwork, v::Vector{Int}) = (h.vec_int4 = v)
cleaned!(h::HybridNetwork, b::Bool) = (h.boolg2 = b)


"""
    Quartet

type that saves the information on a given 4-taxon subset. It contains the following attributes:

- `number`: integer
- `taxon`: vector of taxon names, like t1 t2 t3 t4
- `obsCF`: vector of observed CF, in order 12|34, 13|24, 14|23
- `logPseudoLik`: log pseudolikelihood of the quartet. 0.0 by default
- `ngenes`: number of gene trees used to compute the observed CF; -1.0 if unknown
- `deltaCF`: The sum of absolute differences between observed and expected CFs
- `sampled`: A boolean denoting whether the quartet is used in computing the likelihood
- `uninformative`: A boolean denoting whether the quartet is not sampled due to being uninformative

see also: [`PhyloNetworks.QuartetT`](@extref) for quartet with data of user-defined type `T`,
using a mapping between quartet indices and quartet taxa.
"""
mutable struct Quartet <: AQuartet
    number::Int
    taxon::Array{String,1} # taxa 1234
    obsCF::Array{Float64,1} # three observed CF in order 12|34, 13|24, 14|23
    expCF::Array{Float64,1} # three expected CF in order 12|34, 13|24, 14|23
    logPseudoLik::Float64 # log pseudolik value for the quartet. 0.0 by default
    ngenes::Float64 # number of gene trees used to compute the obsCV, default -1.; Float in case ngenes is average
    deltaCF::Float64 # sum of absolute differences of obsCF - expCF 
    sampled::Bool # false if quartet is not sampled for network optimization, default true
    uninformative::Bool # true if quartet is not sampled because it failed qinfTest, default false
    # inner constructor: to guarantee obsCF are only three and add up to 1
    function Quartet(number::Integer,t1::AbstractString,t2::AbstractString,t3::AbstractString,t4::AbstractString,obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99 < sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        new(number,[t1,t2,t3,t4],obsCF,[],0.0,-1.0, 0.0, true, false);
    end
    function Quartet(number::Integer,t1::Array{String,1},obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99< sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        size(t1,1) != 4 ? error("array of taxa should have size 4, not $(size(t1,1))") : nothing
        0.0 <= obsCF[1] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[1]) for $(t1)")
        0.0 <= obsCF[2] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[2]) for $(t1)")
        0.0 <= obsCF[3] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[3]) for $(t1)")
        new(number,t1,obsCF,[],0.0,-1.0, 0.0, true, false);
    end
    Quartet() = new(0,[],[],[],0.0,-1.0, 0.0, true, false)
end


# Data on quartet concordance factors -------

"""
    DataCF

type that contains the following attributes:

- `quartet` (vector of Quartets)
- `numQuartets`
- `tree` (vector of trees: empty if a table of CF was input instead of list of trees)
- `numTrees` (-1 if a table CF was input instead of list of trees)
- `repSpecies` (taxon names that were repeated in table of CF or input gene trees: used inside snaq for multiple alleles case)

The list of `Quartet` may be accessed with the attribute `.quartet`.
If the input was a list of trees, the `HybridNetwork`'s can be accessed with the attribute `.tree`.
For example, if the `DataCF` object is named `d`, `d.quartet[1]` will show the first quartet
and `d.tree[1]` will print the first input tree.
"""
mutable struct DataCF # fixit
    quartet::Array{Quartet,1} # array of quartets read from CF output table or list of quartets in file
    numQuartets::Integer # number of quartets
    tree::Vector{HybridNetwork} #array of input gene trees
    numTrees::Integer # number of gene trees
    repSpecies::Vector{String} #repeated species in the case of multiple alleles
    DataCF(quartet::Array{Quartet,1}) = new(quartet,length(quartet),[],-1,[])
    DataCF(quartet::Array{Quartet,1},trees::Vector{HybridNetwork}) = new(quartet,length(quartet),trees,length(trees),[])
    DataCF() = new([],0,[],-1,[])
end

# aux type for the updateBL function
mutable struct EdgeParts
    edgenum::Int
    part1::Vector{Node}
    part2::Vector{Node}
    part3::Vector{Node}
    part4::Vector{Node}
end

# Pretty-printing for custom structs
function Base.show(io::IO,d::DataCF)
    print(io,"Object DataCF\n")
    print(io,"number of quartets: $(d.numQuartets)\n")
    if(d.numTrees != -1)
        print(io,"number of trees: $(d.numTrees)\n")
    end
end

function Base.show(io::IO,q::Quartet)
    print(io,"number: $(q.number)\n")
    print(io,"taxon names: $(q.taxon)\n")
    print(io,"observed CF: $(q.obsCF)\n")
    print(io,"expected CF: $(round.(q.expCF, digits=8)) (meaningless before estimation)\n")
    print(io,"pseudo-deviance under last used network: $(q.logPseudoLik) (meaningless before estimation)\n")
    if(q.ngenes != -1)
        print(io,"number of genes used to compute observed CF: $(q.ngenes)\n")
    end
end