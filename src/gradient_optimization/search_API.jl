using PhyloNetworks


function search(N::HybridNetwork, q, hmax::Int; maxeval::Int=500, maxequivPLs::Int=100)

    N = deepcopy(N);
    if N.isrooted
        semidirect_network!(N)
        N.isrooted = false
    end

    logPLs = Array{Float64}(undef, maxeval)
    logPLs[1] = compute_logPL(N, q);
    prop_Ns = [N];
    unchanged_iters = 0

    for j = 2:maxeval
        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = propose_topology(N, hmax)
        @debug "Proposed: $(writenewick(Nprime, round=true))"
        
        # 2. Check for identifiability

        # 3. Optimize branch lengths
        q_eqns, _ = find_quartet_equations(Nprime)
        optimize_bls!(Nprime, q_eqns, q)

        # 4. Compute logPL
        Nprime_logPL = compute_logPL(q_eqns, Nprime.edge, q)
        println("\t\t$(Nprime_logPL)")
        push!(prop_Ns, Nprime)

        # 5. Accept / reject
        if Nprime_logPL == NaN || abs(Nprime_logPL) < 1e-100
            error("Nprime_logPL = $(Nprime_logPL)")
        end

        if Nprime_logPL > logPLs[j-1]
            N = Nprime
            logPLs[j] = Nprime_logPL
            unchanged_iters = 0
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

    return N, logPLs

end



"""
Takes network `N` and modifies it with topological moves to generate a new proposal network.
"""
function propose_topology(N::HybridNetwork, hmax::Int)::HybridNetwork

    if N.numhybrids < hmax && rand() < 0.20
        @debug "MOVE: add_random_hybrid!"
        Nprime = deepcopy(N)
        add_random_hybrid!(Nprime)
        return Nprime
    end

    @debug "MOVE: perform_random_rNNI!"
    Nprime = deepcopy(N)
    perform_random_rNNI!(Nprime)
    return Nprime

end