# Observed quartet concordance factors, read off a set of gene trees.

"""
    observedCF4taxaidx(a, b, c, d, lcadepth, ntaxa, ntrees)

Observed quartet CF (order 12|34, 13|24, 14|23) for the taxa with canonical indices
`a < b < c < d`, read out of [`LazyQuartetCF`](@ref)'s pairwise LCA-depth table instead of
by walking gene trees.

In a binary tree the pair whose LCA is deepest is a cherry, and the cherry determines which
topology the tree displays, so each gene tree costs 6 integer lookups. A depth of `-1` marks
a pair with a taxon missing from that tree; such trees are skipped. If no gene tree holds all
4 taxa the CFs are `NaN`, as elsewhere in the package.
"""
function observedCF4taxaidx(a::Int, b::Int, c::Int, d::Int, lcadepth::Matrix{Int16},
                            ntaxa::Int, ntrees::Int)::NTuple{3,Float64}
    pab = pairindex(a, b, ntaxa); pcd = pairindex(c, d, ntaxa)
    pac = pairindex(a, c, ntaxa); pbd = pairindex(b, d, ntaxa)
    pad = pairindex(a, d, ntaxa); pbc = pairindex(b, c, ntaxa)

    n1 = 0; n2 = 0; n3 = 0
    @inbounds for t = 1:ntrees
        dab = lcadepth[t, pab]
        dab < 0 && continue
        dcd = lcadepth[t, pcd]
        dcd < 0 && continue
        dac = lcadepth[t, pac]; dbd = lcadepth[t, pbd]
        dad = lcadepth[t, pad]; dbc = lcadepth[t, pbc]

        m1 = ifelse(dab > dcd, dab, dcd)
        m2 = ifelse(dac > dbd, dac, dbd)
        m3 = ifelse(dad > dbc, dad, dbc)

        if m1 >= m2 && m1 >= m3
            n1 += 1
        elseif m2 >= m3
            n2 += 1
        else
            n3 += 1
        end
    end

    suma = n1 + n2 + n3
    return (n1 / suma, n2 / suma, n3 / suma)
end


"""
    observedCF4taxa(taxa4, trees, chainmaps) -> (obsCF, ngenes)

Observed quartet CF for the 4 taxa in `taxa4` (order 12|34, 13|24, 14|23), computed by
walking each gene tree's topology (see [`whichquartet`](@ref)).

`chainmaps[k]` maps leaf name => [`ancestorchain`](@ref) for `trees[k]`, as built once by
[`LazyQuartetCF`](@ref)'s constructor, so neither the "does this tree hold all 4 taxa" check
nor the tree walk has to re-scan or re-climb a tree per quartet.

If no gene tree contains all 4 taxa, `ngenes == 0` and `obsCF` is `[NaN, NaN, NaN]`.
"""
function observedCF4taxa(taxa4::Vector{String}, trees::Vector{HybridNetwork}, chainmaps::Vector{Dict{String,Vector{Node}}})::Tuple{Vector{Float64},Float64}
    suma = 0
    sum12 = 0
    sum13 = 0
    sum14 = 0
    for (t, chainmap) in zip(trees, chainmaps)
        if all(tax -> haskey(chainmap, tax), taxa4)
            res = whichquartet(chainmap, taxa4)
            if res == 1
                sum12 += 1
            elseif res == 2
                sum13 += 1
            elseif res == 3
                sum14 += 1
            end
            suma += 1
        end
    end
    return [sum12/suma, sum13/suma, sum14/suma], Float64(suma)
end


"""
Adapted version of `trytreelikequartet` that only returns the quartet type (1 ab|cd, 2 ac|bd, or 3 ad|bc)
of an input that is guaranteed to be a tree, NOT a network. `chainmap` is a
`Dict{String,Vector{Node}}` mapping leaf name => `ancestorchain(leaf)` for the tree that
`taxa4` are being looked up in (see [`observedCF4taxa`](@ref)), precomputed once per
gene tree by [`LazyQuartetCF`](@ref) -- so this never needs to scan or climb the tree.
"""
function whichquartet(chainmap::Dict{String,Vector{Node}}, taxa4::AbstractVector{String})::Int
    chaina = chainmap[taxa4[1]]
    chainb = chainmap[taxa4[2]]
    chainc = chainmap[taxa4[3]]
    chaind = chainmap[taxa4[4]]

    path_ab = pathedges(chaina, chainb)
    path_cd = pathedges(chainc, chaind)
    path_ac = pathedges(chaina, chainc)
    path_bd = pathedges(chainb, chaind)

    if aredisjointedges(path_ab, path_cd)
        return 1
    elseif aredisjointedges(path_ac, path_bd)
        return 2
    else
        return 3
    end
end


"""
    cfmatrixorder(dcf::DataCF) -> (qorder, colmaps)

Correspondence between the quartets in `dcf` and the rows and columns of the CF matrices
used internally by SNaQ (see [`gatherCFmatrix`](@ref) and [`computeexpectedCFmatrix`](@ref)):
row `j` corresponds to quartet `dcf.quartet[qorder[j]]`, and column `k` of row `j`
corresponds to entry `colmaps[j][k]` of that quartet's `obsCF` and `expCF` vectors.

Any code moving CFs between a `DataCF` and these matrices must go through this function,
in both directions, or a CF could be paired with the wrong quartet or the wrong split.
Errors if `dcf` is lazy.
"""
function cfmatrixorder(dcf::DataCF)::Tuple{Vector{Int},Vector{SVector{3,Int}}}
    checknotlazy(dcf, "cfmatrixorder")
    # Helper function for more legible code later
    minmax(i1::Int, i2::Int)::Tuple{Int,Int} = (min(i1, i2), max(i1, i2))

    # This sorting function is what we use to take the set of
    # quartets in `dcf` as they appear and quickly determine
    # the rearrangement that SNaQ's API is expecting
    function labelsorter(a::Vector{String}, b::Vector{String})::Bool
        for j = 4:-1:1
            a[j] < b[j] && return true
            b[j] < a[j] && return false
        end
        return false
    end
    qorder = sortperm(dcf.quartet, lt = (a, b) -> labelsorter(sort(a.taxon), sort(b.taxon)))

    # column of `obsCF` for the split that pairs the taxa at positions i1 and i2 in `q.taxon`
    splitcolumn(i1::Int, i2::Int)::Int =
        minmax(i1, i2) in ((1, 2), (3, 4)) ? 1 :
        minmax(i1, i2) in ((1, 3), (2, 4)) ? 2 : 3

    colmaps = Vector{SVector{3,Int}}(undef, length(qorder))
    for (j, qidx) in enumerate(qorder)
        # taxonperm[k]: position in `q.taxon` of the k-th taxon in sorted order
        taxonperm = sortperm(dcf.quartet[qidx].taxon)
        colmaps[j] = SVector(splitcolumn(taxonperm[1], taxonperm[2]),
                             splitcolumn(taxonperm[1], taxonperm[3]),
                             splitcolumn(taxonperm[1], taxonperm[4]))
    end
    return qorder, colmaps
end


"""
Takes a `DataCF` object `dcf` and returns a `Matrix{Float64}`
corresponding to the observed CF values of each quartet
in `dcf` ordered in the way that `SNaQ` expects internally
(see [`cfmatrixorder`](@ref)).
Errors if `dcf` is lazy: gathering would compute the CFs of every quartet.
"""
function gatherCFmatrix(dcf::DataCF)::Matrix{Float64}
    checknotlazy(dcf, "gatherCFmatrix")
    qorder, colmaps = cfmatrixorder(dcf)
    oCF_matrix = zeros(length(dcf.quartet), 3)
    for (j, qidx) in enumerate(qorder)
        oCF_matrix[j, :] .= dcf.quartet[qidx].obsCF[colmaps[j]]
    end
    return oCF_matrix
end


"""
    gatherngenes(dcf::DataCF)

Vector of the number of genes `q.ngenes` of each quartet in `dcf`, in the order of the rows
of [`gatherCFmatrix`](@ref).
"""
gatherngenes(dcf::DataCF)::Vector{Float64} = [dcf.quartet[i].ngenes for i in cfmatrixorder(dcf)[1]]


"""
    setexpectedCFs!(dcf::DataCF, net::HybridNetwork, ρ=0.0)

Stores in `q.expCF`, for each quartet `q` in `dcf`, the CFs expected from `net` with its
current parameters and inheritance correlation `ρ`. Errors if `dcf` is lazy.
"""
function setexpectedCFs!(dcf::DataCF, net::HybridNetwork, ρ::Real=0.0)::DataCF
    qorder, colmaps = cfmatrixorder(dcf)
    eCFs = computeexpectedCFmatrix(net, ρ)
    size(eCFs, 1) == length(qorder) || error("the network has $(size(eCFs, 1)) quartets " *
        "but the DataCF has $(length(qorder)).")
    for (j, qidx) in enumerate(qorder)
        expCF = zeros(3)
        expCF[colmaps[j]] .= eCFs[j, :]
        dcf.quartet[qidx].expCF = expCF
    end
    return dcf
end
