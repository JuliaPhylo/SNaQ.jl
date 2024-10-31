# circularity: a node has a vector of edges, and an edge has a vector of nodes

# type created from a HybridNetwork only to extract a given quartet
"""
    QuartetNetwork(net::HybridNetwork)

Subtype of `Network` abstract type.
A `QuartetNetwork` object is an internal type used to calculate the
expected CFs of quartets on a given network.
Attributes of the `QuartetNetwork` objects need not be updated at a given time (see below).

The procedure to calculate expected CFs for a given network is as follows:
1. A `QuartetNetwork` object is created for each `Quartet` using
   `extractQuartet!(net,d)` for `net::HybridNetwork` and `d::DataCF`
2. The vector `d.quartet` has all the `Quartet` objects, each with a `QuartetNetwork`
   object (`q.qnet`). Attibutes in `QuartetNetwork` are not updated at this point
3. Attributes in `QuartetNetwork` are partially updated when calculating the
   expected CF (`calculateExpCFAll!`). To calculate the expected CF for this quartet,
   we need to update the attributes: `which`, `typeHyb`, `t1`, `split`, `formula`, `expCF`.
   To do this, we need to modify the `QuartetNetwork` object (i.e. merge edges,...).
   But we do not want to modify it directly because it is connected to the original
   `net` via a map of the edges and nodes, so we use a deep copy:
   `qnet=deepcopy(q.qnet)` and then `calculateExpCFAll!(qnet)`.
   Attributes that are updated on the original `QuartetNetwork` object `q.qnet` are:
    - `q.qnet.hasEdge`: array of booleans of length equal to `net.edge` that shows which identifiable edges and gammas of `net` (`net.ht`) are in `qnet` (and still identifiable). Note that the first elements of the vector correspond to the gammas.
    - `q.qnet.index`: length should match the number of trues in `qnet.hasEdge`. It has the indexes in `qnet.edge` from the edges in `qnet.hasEdge`. Note that the first elements of the vector correspond to the gammas.
    - `q.qnet.edge`: list of edges in `QuartetNetwork`. Note that external edges in `net` are collapsed when they appear in `QuartetNetwork`, so only internal edges map directly to edges in `net`
    - `q.qnet.expCF`: expected CF for this `Quartet`


Why not modify the original `QuartetNetwork`? We wanted to keep the original
`QuartetNetwork` stored in `DataCF` with all the identifiable edges, to be able
to determine if this object had been changed or not after a certain optimization.

The process is:

1. Deep copy of full network to create `q.qnet` for `Quartet q`.
   This `QuartetNetwork` object has only 4 leaves now, but does not have merged edges
   (the identifiable ones) so that we can correspond to the edges in net.
   This `QuartetNetwork` does not have other attributes updated.
2. For the current set of branch lengths and gammas, we can update the attributes
   in `q.qnet` to compute the expected CF. The functions that do this will "destroy"
   the `QuartetNetwork` object by merging edges, removing nodes, etc... So, we do
   this process in `qnet=deepcopy(q.qnet)`, and at the end, only update `q.qnet.expCF`.
3. After we optimize branch lengths in the full network, we want to update the
   branch lengths in `q.qnet`. The edges need to be there (which is why we do
   not want to modify this `QuartetNetwork` object by merging edges), and
   we do not do a deep-copy of the full network again. We only change the values
   of branch lengths and gammas in `q.qnet`, and we can re-calculate the expCF
   by creating a deep copy `qnet=deepcopy(q.qnet)` and run the other functions
   (which merge edges, etc) to get the `expCF`.

Future work: there are definitely more efficient ways to do this (without the deep copies).
In addition, currently edges that are no longer identifiable in `QuartetNetwork`
do not appear in `hasEdge` nor `index`. Need to study this.

```jldoctest
julia> net0 = readTopology("(s17:13.76,(((s3:10.98,(s4:8.99,s5:8.99)I1:1.99)I2:0.47,(((s6:2.31,s7:2.31)I3:4.02,(s8:4.97,#H24:0.0::0.279)I4:1.36)I5:3.64,((s9:8.29,((s10:2.37,s11:2.37)I6:3.02,(s12:2.67,s13:2.67)I7:2.72)I8:2.89)I9:0.21,((s14:2.83,(s15:1.06,s16:1.06)I10:1.78)I11:2.14)#H24:3.52::0.72)I12:1.47)I13:1.48)I14:1.26,(((s18:5.46,s19:5.46)I15:0.59,(s20:4.72,(s21:2.40,s22:2.40)I16:2.32)I17:1.32)I18:2.68,(s23:8.56,(s1:4.64,s2:4.64)I19:3.92)I20:0.16)I21:3.98)I22:1.05);");

julia> net = readTopologyLevel1(writeTopology(net0)) ## need level1 attributes for functions below
HybridNetwork, Un-rooted Network
46 edges
46 nodes: 23 tips, 1 hybrid nodes, 22 internal tree nodes.
tip labels: s17, s3, s4, s5, ...
(s4:8.99,s5:8.99,(s3:10.0,((((s6:2.31,s7:2.31)I3:4.02,(s8:4.97,#H24:0.0::0.279)I4:1.36)I5:3.64,((s9:8.29,((s10:2.37,s11:2.37)I6:3.02,(s12:2.67,s13:2.67)I7:2.72)I8:2.89)I9:0.21,((s14:2.83,(s15:1.06,s16:1.06)I10:1.78)I11:2.14)#H24:3.52::0.721)I12:1.47)I13:1.48,((((s18:5.46,s19:5.46)I15:0.59,(s20:4.72,(s21:2.4,s22:2.4)I16:2.32)I17:1.32)I18:2.68,(s23:8.56,(s1:4.64,s2:4.64)I19:3.92)I20:0.16)I21:3.98,s17:10.0)I22:1.26)I14:0.47)I2:1.99)I1;

julia> q1 = Quartet(1,["s1", "s16", "s18", "s23"],[0.296,0.306,0.398])
number: 1
taxon names: ["s1", "s16", "s18", "s23"]
observed CF: [0.296, 0.306, 0.398]
pseudo-deviance under last used network: 0.0 (meaningless before estimation)
expected CF under last used network: Float64[] (meaningless before estimation)

julia> qnet = SNaQ.extractQuartet!(net,q1)
taxa: ["s1", "s16", "s18", "s23"]
number of hybrid nodes: 1

julia> sum([e.istIdentifiable for e in net.edge]) ## 23 identifiable edges in net
23

julia> idedges = [ee.number for ee in net.edge[[e.istIdentifiable for e in net.edge]]];

julia> print(idedges)
[5, 6, 9, 11, 12, 13, 17, 20, 21, 22, 26, 27, 28, 29, 30, 31, 34, 38, 39, 40, 44, 45, 46]

julia> length(qnet.hasEdge) ## 24 = 1 gamma + 23 identifiable edges
24

julia> sum(qnet.hasEdge) ## 8 = 1 gamma + 7 identifiable edges in qnet
8

julia> print(idedges[qnet.hasEdge[2:end]]) ## 7 id. edges: [12, 13, 29, 30, 31, 45, 46]
[12, 13, 29, 30, 31, 45, 46]

julia> qnet.edge[qnet.index[1]].number ## 11 = minor hybrid edge
11
```
"""
mutable struct QuartetNetwork <: Network
    numTaxa::Int
    numNodes::Int
    numEdges::Int
    node::Array{Node,1}
    edge::Array{Edge,1}
    hybrid::Array{Node,1} # array of hybrid nodes in network
    leaf::Array{Node,1} # array of leaves
    numHybrids::Int # number of hybrid nodes
    hasEdge::Array{Bool,1} # array of boolean with all the original identifiable edges of HybridNetwork and gammas (net.ht)
    quartetTaxon::Array{String,1} # the quartet taxa in the order it represents. Points to same array as its Quartet.taxon
    which::Int8 # 0 it tree quartet, 1 is equivalent to tree quartet and 2 if two minor CF different, default -1
    typeHyb::Array{Int8,1} #array with the type of hybridization of each hybrid node in the quartet
    t1::Float64 # length of internal edge, used when qnet.which=1, default = -1
    names::Array{String,1} # taxon and node names, same order as in network.node
    split::Array{Int8,1} # split that denotes to which side each leaf is from the split, i.e. [1,2,2,1] means that leaf1 and 4 are on the same side of the split, default -1,-1,-1,-1
    formula::Array{Int8,1} # array for qnet.which=1 that indicates if the expCf is major (1) or minor (2) at qnet.expCF[i] depending on qnet.formula[i], default -1,-1,-1
    expCF::Array{Float64,1} # three expected CF in order 12|34, 13|24, 14|23 (matching obsCF from qnet.quartet), default [0,0,0]
    indexht::Vector{Int} # index in net.ht for each edge in qnet.ht
    changed::Bool # true if the expCF would be changed with the current parameters in the optimization, to recalculate, default true
    index::Vector{Int} # index in qnet.edge (qnet.node for gammaz) of the members in qnet.indexht to know how to find quickly in qnet
    # inner constructor
    function QuartetNetwork(net::HybridNetwork)
        net2 = deepcopy(net); #fixit: maybe we dont need deepcopy of all, maybe only arrays
        new(net2.numTaxa,net2.numNodes,net2.numEdges,net2.node,net2.edge,net2.hybrid,net2.leaf,net2.numHybrids, [true for e in net2.edge],[],-1,[], -1.,net2.names,Int8[-1,-1,-1,-1],Int8[-1,-1,-1],[0,0,0],[],true,[])
    end
    QuartetNetwork() = new(0,0,0,[],[],[],[],0,[],[],-1,[],-1.0,[],[],[],[],[],true,[])
end

abstract type AQuartet end

"""
    Quartet

type that saves the information on a given 4-taxon subset. It contains the following attributes:

- number: integer
- taxon: vector of taxon names, like t1 t2 t3 t4
- obsCF: vector of observed CF, in order 12|34, 13|24, 14|23
- logPseudoLik
- ngenes: number of gene trees used to compute the observed CF; -1.0 if unknown
- qnet: [`QuartetNetwork`](@ref), which saves the expCF after snaq estimation to
  emphasize that the expCF depend on a specific network, not the data

see also: [`QuartetT`](@ref) for quartet with data of user-defined type `T`,
using a mapping between quartet indices and quartet taxa.
"""
mutable struct Quartet <: AQuartet
    number::Int
    taxon::Array{String,1} # taxa 1234. qnet.quartetTaxon points to the same array.
    obsCF::Array{Float64,1} # three observed CF in order 12|34, 13|24, 14|23
    qnet::QuartetNetwork # quartet network for the current network (want to keep as if private attribute)
    logPseudoLik::Float64 # log pseudolik value for the quartet. 0.0 by default
    ngenes::Float64 # number of gene trees used to compute the obsCV, default -1.; Float in case ngenes is average
    deltaCF::Float64 # sum of absolute differences of obsCF - expCF 
    sampled::Bool # false if quartet is not sampled for network optimization, default true
    uninformative::Bool # true if quartet is not sampled because it failed qinfTest, default false
    # inner constructor: to guarantee obsCF are only three and add up to 1
    function Quartet(number::Integer,t1::AbstractString,t2::AbstractString,t3::AbstractString,t4::AbstractString,obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99 < sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        new(number,[t1,t2,t3,t4],obsCF,QuartetNetwork(),0.0,-1.0, 0.0, true, false);
    end
    function Quartet(number::Integer,t1::Array{String,1},obsCF::Array{Float64,1})
        size(obsCF,1) != 3 ? error("observed CF vector should have size 3, not $(size(obsCF,1))") : nothing
        0.99< sum(obsCF) < 1.02 || @warn "observed CF should add up to 1, not $(sum(obsCF))"
        size(t1,1) != 4 ? error("array of taxa should have size 4, not $(size(t1,1))") : nothing
        0.0 <= obsCF[1] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[1]) for $(t1)")
        0.0 <= obsCF[2] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[2]) for $(t1)")
        0.0 <= obsCF[3] <= 1.0 || error("obsCF must be between (0,1), but it is $(obsCF[3]) for $(t1)")
        new(number,t1,obsCF,QuartetNetwork(),0.0,-1.0, 0.0, true, false);
    end
    Quartet() = new(0,[],[],QuartetNetwork(),0.0,-1.0, 0.0, true, false)
end

"""
QuartetT{T}

Generic type for 4-taxon sets. Fields:
- `number`: rank of the 4-taxon set
- `taxonnumber`: static vector of 4 integers, assumed to be distinct and sorted
- `data`: object of type `T`

For easier look-up, a unique mapping is used between the rank (`number`) of a
4-taxon set and its 4 taxa (see [`quartetrank`](@ref) and [`nchoose1234`](@ref)):

rank-1 = (t1-1) choose 1 + (t2-1) choose 2 + (t3-1) choose 3 + (t4-1) choose 4

# examples

```jldoctest
julia> nCk = SNaQ.nchoose1234(5)
6×4 Matrix{Int64}:
 0   0   0  0
 1   0   0  0
 2   1   0  0
 3   3   1  0
 4   6   4  1
 5  10  10  5

julia> SNaQ.QuartetT(1,3,4,6, [.92,.04,.04, 100], nCk)
4-taxon set number 8; taxon numbers: 1,3,4,6
data: [0.92, 0.04, 0.04, 100.0]
```
"""
struct QuartetT{T} <: AQuartet where T
    number::Int
    taxonnumber::StaticArrays.SVector{4,Int}
    data::T
end
function Base.show(io::IO, obj::QuartetT{T}) where T
    disp = "4-taxon set number $(obj.number); taxon numbers: "
    disp *= join(obj.taxonnumber,",")
    disp *= "\ndata: "
    print(io, disp)
    print(io, obj.data)
end
function QuartetT(tn1::Int,tn2::Int,tn3::Int,tn4::Int, data::T, nCk::Matrix, checksorted=true::Bool) where T
    if checksorted
        (tn1<tn2 && tn2<tn3 && tn3<tn4) || error("taxon numbers must be sorted")
    end
    QuartetT{T}(quartetrank(tn1,tn2,tn3,tn4,nCk), SVector(tn1,tn2,tn3,tn4), data)
end

"""
    quartetrank(t1,t2,t3,t4, nCk::Matrix)
    quartetrank([t1,t2,t3,t4], nCk)

Return the rank of a four-taxon set with taxon numbers `t1,t2,t3,t4`,
assuming that `ti`s are positive integers such that t1<t2, t2<t3 and t3<t4
(assumptions not checked!).
`nCk` should be a matrix of "n choose k" binomial coefficients:
see [`nchoose1234`](@ref).

# examples

```jldoctest
julia> nCk = SNaQ.nchoose1234(5)
6×4 Matrix{Int64}:
 0   0   0  0
 1   0   0  0
 2   1   0  0
 3   3   1  0
 4   6   4  1
 5  10  10  5

julia> SNaQ.quartetrank([1,2,3,4], nCk)
1

julia> SNaQ.quartetrank([3,4,5,6], nCk)
15
```
"""
@inline function quartetrank(tnum::AbstractVector, nCk::Matrix)
    quartetrank(tnum..., nCk)
end
@inline function quartetrank(t1::Int, t2::Int, t3::Int, t4::Int, nCk::Matrix)
    # rank-1 = t1-1 choose 1 + t2-1 choose 2 + t3-1 choose 3 + t4-1 choose 4
    return nCk[t1,1] + nCk[t2,2] + nCk[t3,3] + nCk[t4,4] + 1
end

"""
    nchoose1234(nmax)

`nmax+1 x 4` matrix containing the binomial coefficient
"n choose k" in row `n+1` and column `k`. In other words,
`M[i,k]` gives "i-1 choose k". It is useful to store these
values and look them up to rank (a large number of) 4-taxon sets:
see [`quartetrank`](@ref).
"""
function nchoose1234(nmax::Int)
    # compute nC1, nC2, nC3, nC4 for n in [0, nmax]: used for ranking quartets
    M = Matrix{Int}(undef, nmax+1, 4)
    for i in 1:(nmax+1)
        M[i,1] = i-1 # n choose 1 = n. row i is for n=i-1
    end
    M[1,2:4] .= 0 # 0 choose 2,3,4 = 0
    for i in 2:(nmax+1)
        for k in 2:4 # to choose k items in 1..n: the largest could be n, else <= n-1
            M[i,k] = M[i-1,k-1] + M[i-1,k]
        end
    end
    return M
end

# Data on quartet concordance factors -------

"""
    DataCF

type that contains the following attributes:

- quartet (vector of Quartets)
- numQuartets
- tree (vector of trees: empty if a table of CF was input instead of list of trees)
- numTrees (-1 if a table CF was input instead of list of trees)
- repSpecies (taxon names that were repeated in table of CF or input gene trees: used inside snaq for multiple alleles case)

The list of Quartet may be accessed with the attribute .quartet.
If the input was a list of trees, the HybridNetwork's can be accessed with the attribute .tree.
For example, if the DataCF object is named d, d.quartet[1] will show the first quartet
and d.tree[1] will print the first input tree.
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

