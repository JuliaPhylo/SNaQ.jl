using PhyloNetworks, SNaQ, DataFrames, Random, Distributed, PhyloCoalSimulations, StatsBase
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("opt_API.jl")
include("search_API.jl")
include("../../test/test_inplace_updates/misc.jl")
include("bfs_search.jl")

using Plots


function BFS_distributed(
    starting_pool::Vector{HybridNetwork},
    d::Int, # depth for each strand to traverse
    l::Int, # maxequivPLs for each strand
    maxfail::Int,   # max number of times a topology can yield no improvement in search before it is discarded
    maxsuccess::Int,# max number of times a top can yield an IMPROVEMENT before it is discarded
    q,              # quartet data
    hmax::Int,
    truenet::HybridNetwork; # JUST HERE FOR TESTING
    lstep::Int=1,   # size by which to reduce the pool when loglik is not improved w/in k
    maxpoolsize::Int=6*starting_pool[1].numtaxa,   # max number of concurrent candidates
    ftolabs::Float64=1e-6,  # optimization param
    ftolrel::Float64=0.05,  # optimization param
    seed::Int=abs(rand(Int) % 10000),
    propQuartets::Float64=1.0,
    maxmultiplicity::Int=maxpoolsize,
    searchargs...   # arguments passed to the search algo
)
    (d == 1 && propQuartets == 1.0) || error("--")

    rng = Random.seed!(seed)
    q_idxs = sample_qindices(starting_pool[1], propQuartets, rng)
    
    # Semi-direct input networks
    for n in starting_pool
        semidirect_network!(n)
    end

    # Setup pool variables
    pool::Vector{HybridNetwork} = Vector{HybridNetwork}([deepcopy_network(n) for n in starting_pool])
    for n in pool
        loglik!(n, optimize_bls!(n, q))
    end
    pool_eqns::Vector{Union{Nothing,Array{QuartetData}}} = [find_quartet_equations(n, q_idxs)[1] for n in pool]
    nfailures::Vector{Int} = fill(0, length(pool))
    nsuccesses::Vector{Int} = fill(0, length(pool))
    multiplicity::Vector{Int} = fill(1, length(pool))
    forceremoved::BitVector = falses(length(pool))
    get_effps() = sum(m for (m, nf, fr) in zip(multiplicity, nfailures, forceremoved) if !fr && m > nf; init=0.0)
    get_vidxs() = findall(nfailures .< maxfail .&& nfailures .< multiplicity .&& .!forceremoved .&& nsuccesses .< maxsuccess)

    tPL = optimize_bls!(truenet, q)
    p = plot(layout=2, size=(800, 400))
    hline!(p[1], [tPL], color=:red, linestyle=:dash, label=false)
    hline!(p[2], [0.0], color=:black, linestyle=:dash, label=false)
    ############# FOR NOW: implement on 1 worker #############
    iter_ii = 0
    while any(j -> nfailures[j] < multiplicity[j] && nfailures[j] < maxfail && !forceremoved[j] && nsuccesses[j] < maxsuccess, 1:length(pool))

        ################################################################################################
        iter_ii += 1
        effpoolsize = get_effps()
        valid_idxs = get_vidxs()
        bestnet = pool[findmax(n -> loglik(n), pool)[2]]
        bestPL = loglik(bestnet)
        worstPL = minimum([loglik(n) for n in pool[valid_idxs]])
        print("\riter_ii=$(iter_ii), effpoolsize=$(effpoolsize), ")
        print("poolsize=$(length(pool)), ")
        print("mult=($(round(mean(multiplicity), digits=2)), $(maximum(multiplicity))), ")
        print("fails=$(sum(nfailures)), ")
        print("worst -loglik=$(round(worstPL, digits=5)), ")
        print("best -loglik=$(round(bestPL, digits=5))      ")


        if iter_ii > 1
            scatter!(p[1], [iter_ii-1], [bestPL], color=:cyan, label=(iter_ii==2 ? "Best PL" : false), markersize=2)
            scatter!(p[1], [iter_ii-1], [worstPL], color=:red, label=(iter_ii==2 ? "Worst PL" : false), markersize=2)
            scatter!(p[2], [iter_ii-1], [hardwiredclusterdistance(truenet, bestnet, false)], color=:green, label=(iter_ii==2 ? "HWCD" : false), markersize=2)
            display(p)
        end
        ################################################################################################


        # Sample from pool
        valid_idxs = get_vidxs()
        sample_idx = sample(rng, valid_idxs, compute_weights(pool[valid_idxs], multiplicity[valid_idxs]))
        net = pool[sample_idx]

        # Perform single iter of search
        next_net, next_eqns = BFS_single_iter(
            rng, net, pool_eqns[sample_idx], d, q, q_idxs, hmax, ftolabs, ftolrel; maxequivPLs=l, searchargs...
        )
        
        # If improvement found: add it to the pool
        # If NO improvement found: increment `nfailures` 
        if next_net !== nothing
            nsuccesses[sample_idx] += 1
            if nsuccesses[sample_idx] >= maxsuccess
                forceremoved[sample_idx] = true
                pool_eqns[sample_idx] = nothing
            end

            match_idx = findfirst(j -> pool[j] ≊ next_net, 1:length(pool))
            if match_idx === nothing
                push!(pool, next_net)
                push!(pool_eqns, next_eqns)
                push!(nfailures, 0)
                push!(nsuccesses, 0)
                push!(multiplicity, 1)
                push!(forceremoved, false)
            else
                # Finding a network that has already been discovered is considered a failure!
                nfailures[sample_idx] += 1
                multiplicity[match_idx] = min(multiplicity[match_idx] + 1, maxmultiplicity)
            end
        else
            nfailures[sample_idx] += 1
        end
        if nfailures[sample_idx] >= maxfail
            forceremoved[sample_idx] = true
            pool_eqns[sample_idx] = nothing
        end

        # Trim pool back down to `maxpoolsize`
        effpoolsize = get_effps()
        if effpoolsize > maxpoolsize
            valid_idxs = get_vidxs()
            while effpoolsize > maxpoolsize
                # find worst entry in pool
                worst_vidx = -1;
                worst_idx = -1;
                worst_PL = Inf;
                for (j, vidx) in enumerate(valid_idxs)
                    if loglik(pool[vidx]) < worst_PL
                        worst_PL = loglik(pool[vidx])
                        worst_vidx = vidx
                        worst_idx = j
                    end
                end
                forceremoved[worst_vidx] = true
                deleteat!(valid_idxs, worst_idx)
                effpoolsize -= multiplicity[worst_vidx]
            end
        end

        # Remove any entries with a less than 1 / (100 * maxpoolsize) chance of being selected
        valid_idxs = get_vidxs()
        length(valid_idxs) == 0 && break
        sample_idx = sample(rng, valid_idxs, compute_weights(pool[valid_idxs], multiplicity[valid_idxs]))
        W = compute_weights(pool[valid_idxs], multiplicity[valid_idxs])
    end
    ##########################################################
    
    best_idx = findmax(n -> loglik(n), pool)[2]
    return pool[best_idx]
end


function BFS_single_iter(
    rng::TaskLocalRNG,
    N::HybridNetwork,
    net_eqns::Array{QuartetData},
    d::Int,
    q,
    q_idxs::AbstractVector{Int},
    hmax::Int,
    ftolabs::Float64,
    ftolrel::Float64;
    restrictions::Function=no_restrictions(),
    α::Real=Inf,
    propQuartets::Real=1.0,
    maxequivPLs::Int=100,
    opt_maxeval::Int=10
)

    d == 1 || error("ONLY IMPLEMENTED WITH d=1 RIGHT NOW")
    propQuartets == 1.0 || error("PROPQUARTETS != 1.0 NOT IMPLEMENTED YET")
    # Parameter bounds enforcement
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    1 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")

    # Copy input net
    N = deepcopy_network(N);
    semidirect_network!(N);

    # Perform the search
    moves_attempted::Vector{<:Tuple} = Vector{Tuple}();   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
                            # here so that we don't redo moves
    for j = 1:maxequivPLs
        # 1. Propose a new topology
        Nprime, prop_move, prop_params = propose_topology(rng, N, hmax, moves_attempted, restrictions);

        # 2. Optimize that topology
        Nprime_logPL, Nprime_eqns = optimize_topology!(
            Nprime, net_eqns, prop_move, prop_params, q, q_idxs,
            opt_maxeval, N.numhybrids != Nprime.numhybrids, rng, α
        )
        loglik!(Nprime, Nprime_logPL)

        # 3. Accept/reject the topology
        loglik(Nprime) - loglik(N) > ftolabs || continue

        # 4. We accepted the topology, so we need to remove bad hybrids if they exist
        bad_Hs = []
        for H in Nprime.hybrid
            if getparentedgeminor(H).gamma <= 0.001
                bad_H = H
            end
        end
        if length(bad_Hs) != 0
            for H in bad_Hs
                remove_hybrid!(Nprime, H)
            end
            Nprime_eqns, _ = find_quartet_equations(Nprime, q_idxs)
        end
        return Nprime, Nprime_eqns
    end

    return nothing, nothing

end


function propose_topology(rng::TaskLocalRNG, N::HybridNetwork, hmax::Int, moves_attempted::Vector{<:Tuple}, restrictions::Function)
    j = 0
    while true
        j += 1
        j >= 1000 && return nothing, nothing, nothing   # couldn't find a valid move
    
        # Copy network
        Nprime = readnewick(writenewick(N));

        # Find a valid move and perform the move
        prop_move, prop_params = generate_move_proposal(Nprime, moves_attempted, hmax, rng)
        apply_move!(Nprime, prop_move, prop_params)
        push!(moves_attempted, (prop_move, prop_params))

        # Check that this move meet identifiability requirements
        removedegree2nodes!(Nprime)
        while shrink2cycles!(Nprime) || shrink3cycles!(Nprime) continue end
        semidirect_network!(Nprime)

        # Check if network meets restrictions
        if !restrictions(Nprime) continue end
        return Nprime, prop_move, prop_params   # successfully found a proposal network
    end
    return nothing, nothing, nothing    # this should be unreachable, but just in case
end


function compute_weights(N::Vector{HybridNetwork}, multiplicity::Vector{Int})::Weights{Float64}
    L::Vector{Float64} = [loglik(n) for n in N]
    sumL = sum(L)
    return Weights([m * max(sumL / l, 0.0) for (l, m) in zip(L, multiplicity)])  # weight: (xi / sum(xi)) ^ {-1}
                                            # b/c smaller magnitude is better
end


function ≊(N1::HybridNetwork, N2::HybridNetwork)::Bool
    if loglik(N1) ≈ 0.0
        loglik(N2) ≈ 0.0 && return true
    end
    abs((loglik(N1) - loglik(N2)) / loglik(N1)) ≤ 0.01 || return false
    return hardwiredclusterdistance(N1, N2, false) == 0
end




net = generate_net(10, 2, 44);
while shrink2cycles!(net) || shrink3cycles!(net) continue end
net.numhybrids, istreechild(net), getlevel(net)
gts = simulatecoalescent(net, 10000, 1);
q, t = countquartetsintrees(gts; showprogressbar=false);
semidirect_network!(net);

rt = @elapsed results = BFS_distributed(gts[1:10], 1, 25, 1, 100, q, net.numhybrids, net; restrictions=restriction_set(; require_strongly_tree_child=true))

istreechild(results)[3] == true
