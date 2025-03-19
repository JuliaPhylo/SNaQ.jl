using PhyloNetworks, Random, Distributed
include("CF_equations.jl")
include("opt_API.jl")
include("../network_moves/add_remove_retic.jl")
include("../network_moves/rNNI_moves.jl")
include("../network_moves/rSPR_moves.jl")
include("../network_moves/move_origin_target.jl")
include("../network_moves/identifiability_moves.jl")
include("../network_properties/network_properties.jl")


API_WARNED = false
"""
API NOT FINALIZED
"""
function multi_search(N::HybridNetwork, q, hmax::Int; runs::Int=10, seed::Int=abs(rand(Int) % 100000), kwargs...)
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")

    global API_WARNED
    if !API_WARNED
        @warn "API NOT YET FINALIZED"
        API_WARNED = true
    end

    # Prep data
    N = readnewick(writenewick(N)); N.isrooted = false;

    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    # Do the runs distributed
    nets_and_PLs = Distributed.pmap(1:runs) do j
        return search(N, q, hmax; seed = run_seeds[j], kwargs...)
    end

    # Consolidate return data
    all_nets = Array{HybridNetwork}(undef, runs)
    all_logPLs = zeros(Float64, runs)
    for j = 1:runs
        all_nets[j], all_logPLs[j] = nets_and_PLs[j]
    end

    # Return
    sort_idx = sortperm(all_logPLs, rev=true)
    return all_nets[sort_idx[1]], all_nets[sort_idx], all_logPLs[sort_idx]
end


"""
API NOT FINALIZED
"""
function search(
    N::HybridNetwork, q, hmax::Int;
    restrictions::Function=no_restrictions(),
    α::Real=Inf,
    propQuartets::Real=1.0,
    maxeval::Int=1500,
    maxequivPLs::Int=200,
    opt_maxeval::Int=25,
    seed::Int=abs(rand(Int) % 100000)
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    1 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")

    # Set the seed
    rng = Random.seed!(seed)

    # TODO: remove this once API finalized
    global API_WARNED
    if !API_WARNED
        @warn "API NOT YET FINALIZED"
        API_WARNED = true
    end

    N = readnewick(writenewick(N));
    semidirect_network!(N)
    perform_random_rNNI!(N, rng)

    logPLs = Array{Float64}(undef, maxeval)
    N_qdata, _, N_params, _ = find_quartet_equations(N);
    logPLs[1] = compute_loss(N_qdata, N_params, q, α);
    unchanged_iters = 0

    moves_proposed = zeros(9)
    moves_accepted = zeros(9)

    for j = 2:maxeval
        # print("\rCurrent best -logPL: $(-round(logPLs[j-1], digits=2))          ($(j)/$(maxeval))                 ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));
        move_proposed = propose_topology!(Nprime, hmax, rng)
        moves_proposed[move_proposed] += 1
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        # 2. Check for identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        # shrink3cycles!(Nprime);
        shrink2cycles!(Nprime);
        # shrink_bad_diamonds!(Nprime);

        # 3. Immediately throw away networks that don't meet restrictions
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            logPLs[j] = logPLs[j-1]
            continue
        end

        # 3. Optimize branch lengths
        @debug "\tgathering quartets"
        Nprime_qdata, _ = find_quartet_equations(Nprime);
        @debug "\toptimizing BLs"
        Nprime_logPL = optimize_bls!(Nprime, Nprime_qdata, q; maxeval=opt_maxeval)
        @debug "Optimized network: $(writenewick(Nprime, round=true))"

        # 4. Remove hybrids with γ ≈ 0
        bad_H = nothing
        for H in N.hybrid
            if getparentedgeminor(H).gamma <= 0.001
                bad_H = H
                break
            end
        end
        if bad_H !== nothing
            @debug "Found H with γ=$(getparentedge(bad_H).gamma), removing it."
            remove_hybrid!(bad_H, N)
        end

        # 5. Accept / reject
        if Nprime_logPL === NaN || abs(Nprime_logPL) < 1e-100
            error("Nprime_logPL = $(Nprime_logPL)")
        end

        if Nprime_logPL - logPLs[j-1] > 1e-8
            N = Nprime
            logPLs[j] = Nprime_logPL
            unchanged_iters = 0
            moves_accepted[move_proposed] += 1
        else
            logPLs[j] = logPLs[j-1]
            unchanged_iters += 1
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            logPLs = logPLs[1:j]
            break
        end
    end
    # println("\rBest -logPL discovered: $(-round(logPLs[length(logPLs)], digits=2))")

    return N, logPLs[length(logPLs)]

end


"""
Takes network `N` and modifies it with topological moves to generate a new proposal network.
Also updates `Nprime_eqns` in-place to hold the equations of the newly proposed network.
"""
function propose_topology!(N::HybridNetwork, hmax::Int, rng::TaskLocalRNG)

    if N.numhybrids < hmax && rand(rng) < 0.50
        @debug "MOVE: add_random_hybrid!"
        add_random_hybrid!(N, rng)
        return 2
    end

    if rand(rng) < 0.5

        @debug "MOVE: perform_random_rNNI!"
        return 2 + perform_random_rNNI!(N, rng)
    
    else

        try
            if rand(rng) < 0.5
                @debug "move_random_reticulate_origin!"
                # move_random_reticulate_origin!(N)
                move_random_reticulate_origin_local!(sample(rng, N.hybrid), N, rng)
                return 7
            else
                @debug "move_random_reticulate_target!"
                # move_random_reticulate_target!(N)
                move_random_reticulate_target_local!(sample(rng, N.hybrid), N, rng)
                return 8
            end
        catch
            # Only catches if there were 0 valid moves for whichever above was chosen
        end
    end

    @debug "MOVE: perform_random_rNNI!"
    return 2 + perform_random_rNNI!(N, rng)

end


