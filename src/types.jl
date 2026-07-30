"""
    SNaQscore(network::HybridNetwork)

Gets the cached composite SNaQ score of the provided network. This value is *proportional* to
the *composite deviance* of the network. To compute the composite log-likelihood of the
network, see [`compositeloglik`](@ref).

If [`computeSNaQscore!`](@ref) or [`fitnumericalparameters!`](@ref) have not been run on the
provided network, and the provided network is not the output of [`snaq!`](@ref),
[`readsnaqnetwork`](@ref), or [`readallsnaqnetworks`](@ref), then the return value
will be meaningless.
"""
SNaQscore(h::HybridNetwork) = h.fscore


"""
    SNaQscore!(network, value)

Sets the cached composite SNaQ score of network `h` to `f`.
Higher values indicate better fit.
"""
SNaQscore!(h::HybridNetwork, f::Real) = (h.fscore = f)


"""
    Quartet

type that saves the information on a given 4-taxon subset. It contains the following attributes:

- `number`: integer
- `taxon`: vector of taxon names, like t1 t2 t3 t4
- `obsCF`: vector of observed CF, in order 12|34, 13|24, 14|23
- `expCF`: vector of CF expected from the current network, in the same order.
  warning: this vector may be obsolete
  (see [`computeSNaQscore!`](@ref) or [`fitnumericalparameters!`](@ref))
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