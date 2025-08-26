"""
    multi_search(
        N::HybridNetwork,
        q::Union{DataCF, AbstractArray{Float64}},
        hmax::Int;
        runs::Int=10,
        seed::Int=42,
        kwargs...
    ) -> Tuple{HybridNetwork, Vector{HybridNetwork}, Vector{Float64}}

Performs multiple independent searches to find the best network topology and parameters
that fit the observed quartet concordance factors.

# Arguments
- `N::Union{HybridNetwork, Vector{HybridNetwork}}`: The starting network topology or vector
    of topologies. If a vector, must be length 1 or of length `runs`.
- `q::Union{DataCF, AbstractArray{Float64}}`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `runs::Int=10`: Number of independent search runs.
- `seed::Int=42`: Random seed for reproducibility.
- `kwargs...`: Additional keyword arguments passed to the [`search`](@ref) function.

# Returns
- `best_network::HybridNetwork`: The network with the best (lowest) negative log pseudo-likelihood.
- `all_networks::Vector{HybridNetwork}`: All networks from the runs, sorted by score.
- `all_scores::Vector{Float64}`: Negative log pseudo-likelihood scores for each network.

# Notes
- This function uses distributed computing to perform searches in parallel.
- It returns the best network found across all runs.
"""
function multi_search(
    N::Union{HybridNetwork, AbstractVector{HybridNetwork}},
    q::Union{DataCF, AbstractMatrix{Float64}},
    hmax::Int;
    # Basic arguments
    runs::Int=10,
    seed::Int=42,
    logprefix::String="",
    outgroup::String="none",
    restrictions::Function=defaultrestrictions(),
    kwargs...
)
    # Verify input parameters
    runs > 0 || error("runs must be > 0 (runs = $(runs)).")
    typeof(N) <: Vector{HybridNetwork} && length(N) != 1 && length(N) != runs && error("If N is a vector, it must be length 1 or length equal to runs (length(N) = $(length(N)), runs = $(runs))")
    N = deepcopy(N)

    if typeof(q) <: DataCF
        # If input is DataCF, make sure there's not a name mismatch
        dcf_names::Vector{String} = []
        for quartet in q.quartet
            for tax in quartet.taxon
                if !(tax ∈ dcf_names)
                    push!(dcf_names, tax)
                end
            end
        end

        # If any taxon are in the input but not the inputs but not the quartets, remove them from the inputs
        # If any taxon are NOT in the input but are in the quartets, error
        for T in (typeof(N) <: HybridNetwork ? [N] : N)
            any(tax -> tax ∉ tiplabels(T), dcf_names) && throw(ErrorException("Taxa in DataCF does not match taxa in input."))
            removetaxa = []
            for taxon in tiplabels(T)
                if taxon ∉ dcf_names
                    push!(removetaxa, taxon)
                end
            end
            if length(removetaxa) > 0
                @warn "The following are in the inputs but not the DataCF, these taxa will be deleted: $(removetaxa)"
                for taxon in removetaxa
                    PhyloNetworks.deleteleaf!(T, taxon)
                end
            end
        end

        if typeof(N) <: HybridNetwork
            length(symdiff(tiplabels(N), dcf_names)) == 0 || throw(ErrorException("Taxa in DataCF does not match taxa in input."))
        else
            all(t -> length(symdiff(tiplabels(t), dcf_names)) == 0, N) || throw(ErrorException("Taxa in DataCF does not match taxa in one of the inputs."))
        end
    else
        # If input is AbstractArray{Float64}, make sure its size is (ntaxa choose 4, 3)
        nrow::Int = binomial(typeof(N) <: HybridNetwork ? N.numtaxa : N[1].numtaxa, 4)
        size(q) == (nrow, 3) || error("Input CF matrix should have size ($(nrow), 3), has size $(size(q)) instead.")
    end

    # Verify the starting network inputs
    Ns::Vector{HybridNetwork} = verifystartingtopologies!(N, outgroup, restrictions)

    # Convert q to a Matrix if it is a DataCF
    if typeof(q) <: DataCF
        qtemp::Matrix{Float64} = Matrix{Float64}(undef, q.numQuartets, 3)
        for j = 1:q.numQuartets
            for k = 1:3
                qtemp[j, k] = q.quartet[j].obsCF[k]
            end
        end
        q = qtemp
        qtemp = Matrix{Float64}(undef, 0, 0)
    end

    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    # Do the runs distributed
    nets_and_PLs = pmap(
        j -> search(
            length(Ns) == 1 ? Ns[1] : Ns[j],
            q, hmax; seed = run_seeds[j], restrictions=restrictions,
            logfile = logprefix == "" ? "" : "$(logprefix)$(run_seeds[j])",
            outgroup=outgroup, kwargs...
        ),
        1:runs
    )

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
Verifies that starting topology(ies) `N` are ready to be optimized.
Modifies the network(s) `N` in-place.
"""
function verifystartingtopologies!(N::Union{HybridNetwork, AbstractVector{HybridNetwork}}, outgroup::String, restrictions::Function)::Vector{HybridNetwork}
    # Copy the input networks
    Ns::Vector{HybridNetwork} = typeof(N) <: HybridNetwork ? [deepcopy_network(N)] : [deepcopy_network(n) for n in N]
    for (j, n) in enumerate(Ns)
        # Prep data
        semidirect_network!(Ns[j]);

        # Make sure all leaf edges have some length so that code later doesn't error
        for E in Ns[j].edge
            E.length = E.length == -1.0 ? 0.0 : E.length
        end
        for H in Ns[j].hybrid
            if 1 ≥ getparentedge(H).gamma ≥ 0 && 1 ≥ getparentedgeminor(H).gamma ≥ 0 && getparentedge(H).gamma + getparentedgeminor(H).gamma ≈ 1
                continue
            end
            getparentedge(H).gamma = 0.5
            getparentedgeminor(H).gamma = 0.5
        end

        # Make sure starting network meets restrictions if any are provided
        restrictions(Ns[j]) || throw(ArgumentError("Starting topology #$(j) does not meet provided restrictions."))

        # If no outgroup exists, go next
        if outgroup == "none" continue end
        
        # Make sure the outgroup exists in this network
        if !any(L -> L.name == outgroup, Ns[j].leaf)
            throw(ArgumentError("Starting topology #$(j) does not contain the supplied outgroup ($(outgroup))."))
        end

        # Try rooting at outgroup if there is one
        try
            PN.rootatnode!(Ns[j], outgroup)
        catch e
            if typeof(e) <: PN.RootMismatch
                throw(ArgumentError("Starting topology #$(j) contains the outgroup but cannot be rooted at the outgroup."))
            else
                rethrow(e)
            end
        end
    end
    return Ns
end


"""
Helper function to log the message `msg` to file `logfile`.
"""
function log_text(logfile::String, msg::String)
    logfile == "" && return
    # get the current time and format it
    timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
    # open the file in append mode and write the log line
    open(logfile, "a") do io
        println(io, "[$timestamp] $msg")
    end
end


"""
Helper function to log proposed and accepted moves to `logfile`.
"""
function log_moves(logfile::String, moves_prop::Dict, moves_acc::Dict, moves_PL::Dict)
    all_keys = sort(collect(keys(moves_prop)))
    min_width::Int = maximum([max(length(string(k)), length(string(moves_prop[k])), length(string(moves_acc[k]))) for k in all_keys])+2
    function expand(str::String)::String
        ret::String = str
        for j = 1:(min_width - length(str))
            ret *= " "
        end
        return ret
    end

    msg::String = repeat("-", 12) * "MOVE ACCEPTANCE RATES" * repeat("-", 12) * "\n"
    msg *= expand("\tmove:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(k))
    end
    msg *= " |\n" * expand("\tproposals:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(moves_prop[k]))
    end
    msg *= " |\n" * expand("\taccepted:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(moves_acc[k]))
    end
    msg *= " |\n" * expand("\taccept %:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(100 * moves_acc[k] / moves_prop[k], digits=0)))
    end
    msg *= " |\n" * expand("\tmean SUCC ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(mean(moves_PL[k]), digits=2)))
    end
    msg *= " |\n" * expand("\tmax SUCC ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(maximum(moves_PL[k], init=0), digits=2)))
    end
    msg *= " |\n" * expand("\tmean ALL ΔPL:")
    for k in all_keys
        msg *= "| "
        msg *= expand(string(round(sum(moves_PL[k]) / moves_prop[k], digits=2)))
    end
    msg *= " |\n\n"

    log_text(logfile, msg)
end



"""
    search(
        N::HybridNetwork,
        q,
        hmax::Int;
        restrictions::Function=defaultrestrictions(),
        α::Real=Inf,
        propQuartets::Real=1.0,
        preopt::Bool=false,
        prehybprob::Real=0.0,
        prehybattempts::Int=5,
        probST::Real=0.3,
        maxeval::Int=Int(1e8),
        maxequivPLs::Int=1500,
        opt_maxeval::Int=10,
        seed::Int=abs(rand(Int) % 100000),
        verbose::Bool=false,
        logfile::String=""
    ) -> Tuple{HybridNetwork, Float64}

Performs a single search for the optimal network topology with gradient-based optimization
of branch lengths and inheritance probabilities.

# Arguments
- `N::HybridNetwork`: The starting network topology.
- `q`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `restrictions::Function=defaultrestrictions()`: Function to enforce restrictions on the proposed networks.
- `α::Real=Inf`: Dirichlet parameter for gene tree heterogeneity model.
- `propQuartets::Real=1.0`: Proportion of quartets to use during optimization.
- `preopt::Bool=false`: Whether to perform a pre-optimization step.
- `prehybprob::Real=0.0`: Probability of attempting pre-optimization hybrid attachments.
- `prehybattempts::Int=5`: Number of attempts for pre-optimization hybrid attachments.
- `probST::Real=0.3`: Probability of performing a subtree move before searching.
- `maxeval::Int=Int(1e8)`: Maximum number of evaluations.
- `maxequivPLs::Int=1500`: Maximum number of equivalent pseudo-likelihood scores to consider.
- `opt_maxeval::Int=10`: Maximum evaluations for optimization.
- `seed::Int=abs(rand(Int) % 100000)`: Random seed for reproducibility.
- `verbose::Bool=false`: Whether to print verbose output.
- `logfile::String=""`: File to log progress.

# Returns
- `best_network::HybridNetwork`: The network with the best (lowest) negative log pseudo-likelihood.
- `best_score::Float64`: The negative log pseudo-likelihood of the best network.
"""
function search(
    N::HybridNetwork,
    q::Union{DataCF, Matrix{Float64}},
    hmax::Int;
    restrictions::Function=defaultrestrictions(),
    α::Real=Inf,
    propQuartets::Real=1.0,
    #preopt::Bool=true,
    preopt::Bool=false,
    #prehybprob::Real=0.3,
    prehybprob::Real=0.0,
    prehybattempts::Int=5,
    probST::Real=0.3,
    maxeval::Int=Int(1e8),
    maxequivPLs::Int=1500,
    liktolAbs::Float64=1e-8,
    liktolRel::Float64=1e-4,
    opt_maxeval::Int=10,
    seed::Int=abs(rand(Int) % 100000),
    verbose::Bool=false,
    logfile::String="",
    outgroup::String="none",
    optargs...
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    1 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 ≤ probST ≤ 1 || error("probST must be in range [0, 1] (probST = $(probST))")
    0 ≤ prehybprob ≤ 1 || error("prehybprob must be in range [0, 1] (prehybprob = $(prehybprob))")
    0 ≤ prehybattempts ≤ Inf || error("prehybattempts must be ≥ 0 (prehybattempts = $(prehybattempts))")
    outgroup == "none" || any(l -> l.name == outgroup, N.leaf) || error("No taxa in N have taxa name $(outgroup) (outgroup name)")

    # Set the seed
    rng = Random.seed!(seed)

    N = readnewick(writenewick(N));
    semidirect_network!(N)
    restrictions(N) || error("N does not meet restrictions IMMEDIATELY")

    if rand(rng) < probST
        try
            perform_rNNI1!(N, sample_rNNI_parameters(N, 1, rng)...);
        catch
        end
    end
    restrictions(N) || error("N does not meet restrictions after probST")

    # Initial pre-opt search
    if preopt
        @debug "Pre-optimizing"
        N, _ = search(N, q, hmax;
            restrictions=restrictions,
            preopt=false,
            probST=0.0,
            maxeval=1000,
            opt_maxeval=opt_maxeval,
            prehybattempts=0,
            maxequivPLs=50
        )
        restrictions(N) || error("N does not meet restrictions after preopt")
    end

    # Try many different ways of adding hybrids and pick the best
    if prehybattempts > 0 && rand(rng) < prehybprob
        @debug "Attempting pre-opt hybrid attachments"
        attempt_prehybs(N, q, hmax, restrictions, prehybattempts, rng)
        restrictions(N) || error("N does not meet restrictions after prehyb")
    end

    # Data used throughout the optimization process
    q_idxs = sample_qindices(N, propQuartets, rng)
    logPLs::Array{Float64} = Array{Float64}(undef, maxeval)
    neq = find_quartet_equations(N, q_idxs);
    N_eqns::Vector{QuartetData} = neq[1];
    N_numparams::Int = length(neq[3])
    logPLs[1] = optimize_bls!(N, N_eqns, q[q_idxs,:], α; maxeval=opt_maxeval, optargs...)
    unchanged_iters = 0

    moves_attempted = [];   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
    moves_proposed = Dict{Symbol,Int}()
    moves_accepted = Dict{Symbol,Int}()
    moves_logPL = Dict{Symbol,Vector{Float64}}()
    last_move = :nothing

    log_text(logfile, "Entering main loop with -logPL = $(logPLs[1])")
    for j = 2:maxeval
        if j % 100 == 0
            log_moves(logfile, moves_proposed, moves_accepted, moves_logPL)
        end

        verbose && print("\rIteration $(j)/$(maxeval) - in a row=$(unchanged_iters)/$(maxequivPLs)              ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        Nprime = readnewick(writenewick(N));

        prop_move, prop_params = generate_move_proposal(Nprime, moves_attempted, hmax, rng)
        last_move = prop_move
        apply_move!(Nprime, prop_move, prop_params)
        push!(moves_attempted, (prop_move, prop_params))
        @debug "Proposed: $(writenewick(Nprime, round=true))"

        if !haskey(moves_proposed, prop_move) moves_proposed[prop_move] = 0 end
        if !haskey(moves_accepted, prop_move) moves_accepted[prop_move] = 0 end
        if !haskey(moves_logPL, prop_move) moves_logPL[prop_move] = Vector{Float64}([]) end
        moves_proposed[prop_move] += 1

        # 2. Check for identifiability
        @debug "Proposed network level: $(getlevel(Nprime))"
        removedegree2nodes!(Nprime);
        while shrink3cycles!(Nprime) continue end
        while shrink2cycles!(Nprime) continue end   # keep shrinking until there is nothing to shrink

        # 2.2 Try re-rooting at the outgroup - if we can't, throw the network away
        if outgroup != "none"
            try
                rootatnode!(Nprime, outgroup)
            catch e
                if typeof(e) <: PN.RootMismatch
                    @debug "Nprime cannot be rooted at outgroup - skipping."
                    log_text(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (cannot reroot at outgroup)")
                    logPLs[j] = logPLs[j-1]
                    continue
                else
                    rethrow(e)
                end
            end
        end

        # 2.3 After removing some edges above, the root may have 2 edge now instead of 3 - we fix that here
        semidirect_network!(Nprime)

        # 3. Immediately throw away networks that don't meet restrictions 
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            log_text(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (restrictions not met)")
            logPLs[j] = logPLs[j-1]
            continue
        end

        # Check whether we can do in-place updates here.
        # We CANNOT do inplace updates if:
        # 1. the number of hybrids changes, OR
        # 2. the number of optimization parameters in the network changed
        Nprime_np::Int = gather_optimization_info(Nprime)[1]
        N_np::Int = gather_optimization_info(N)[1]
        cannot_do_inplace::Bool = N.numhybrids != Nprime.numhybrids || N_np != Nprime_np

        # 4. Optimize branch lengths and compute logPL
        Nprime_logPL, Nprime_eqns = optimize_topology!(
            Nprime, N_eqns, prop_move, prop_params, q, q_idxs,
            opt_maxeval, cannot_do_inplace, rng, α; optargs...
        )
        Nprime_logPL == -Inf && error("Nprime_logPL is -Inf?? newick: $(writenewick(Nprime, round=true))\nold network: $(writenewick(N, round=true))\nprop move: $(prop_move)\nprop params: $(prop_params)")
        # compute_loss(Nprime, q) == Nprime_logPL || error("LOGPLS NOT EQUAL AFTER MOVE $(prop_move)")

        # 5. Accept / reject
        (Nprime_logPL === NaN || abs(Nprime_logPL) < eps()) && error("""
            Nprime_logPL = $(Nprime_logPL)
            $(writenewick(Nprime, round=true))
        """)
        if Nprime_logPL - logPLs[j-1] > liktolAbs && (logPLs[j-1] - Nprime_logPL) / logPLs[j-1] > liktolRel
            # Update current topology info
            N = Nprime
            N_eqns = Nprime_eqns
            logPLs[j] = Nprime_logPL
            moves_accepted[prop_move] += 1
            push!(moves_logPL[prop_move], logPLs[j] - logPLs[j-1])

            # Log acceptance
            log_text(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) ACCEPTED $(prop_move), new -logPL=$(round(logPLs[j], digits=6))")

            # Update tracking vars
            unchanged_iters = 0
            moves_attempted = []
        else
            logPLs[j] = logPLs[j-1]
            unchanged_iters += 1

            # Log rejection and reason
            log_text(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) ($(round(Nprime_logPL, digits=3)) < $(round(logPLs[j], digits=3)))")
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            logPLs = logPLs[1:j]
            break
        end
    end

    log_text(logfile, "Search complete.\n\n")
    log_moves(logfile, moves_proposed, moves_accepted, moves_logPL)

    return N, logPLs[length(logPLs)]

end


function attempt_prehybs(net::HybridNetwork, q, hmax::Int, restrictions::Function, attempts::Int, rng::TaskLocalRNG; optargs...)
    best_net = deepcopy_network(net)
    last_net = best_net
    best_PL = optimize_bls!(net, q; optargs...)

    for nhyb = net.numhybrids:hmax
        best_PL = -Inf

        for j = 1:attempts
            local N0::HybridNetwork
            found_valid_move = false
            
            # Try to add a hybrid that complies w/ the given restrictions 100 times
            # If none are found, we just skip
            for restr_attempt = 1:100
                N0 = deepcopy_network(last_net)
                apply_move!(N0, :add_hybrid, Tuple(sample_add_hybrid_parameters(N0, rng)))
                if restrictions(N0)
                    found_valid_move = true
                    break
                end
            end
            found_valid_move || continue

            new_PL = optimize_bls!(N0, q; optargs...)
            if new_PL > best_PL
                best_PL = new_PL
                best_net = N0
            end
        end

        best_net == last_net && break   # Couldn't find a move that is valid under `restrictions`
        last_net = best_net
    end
    return best_net, best_PL
end


"""
Applies the move `move` on parameters `params` to network `N`.
"""
function apply_move!(N::HybridNetwork, move::Symbol, params::Tuple)
    if move == :add_hybrid
        return add_hybrid!(N, params...)
    elseif move == :rNNI1
        return perform_rNNI1!(N, params...)
    elseif move == :rNNI2
        return perform_rNNI2!(N, params...)
    elseif move == :rNNI3
        return perform_rNNI3!(N, params...)
    elseif move == :rNNI4
        return perform_rNNI4!(N, params...)
    elseif move == :retic_origin || move == :retic_origin_local
        return move_reticulate_origin!(N, params...)
    elseif move == :retic_target || move == :retic_target_local
        return move_reticulate_target!(N, params...)
    elseif move == :rSPR
        return perform_rSPR!(N, params...)
    elseif move == :flip_hybrid
        return fliphybrid!(N, params[1])
    end

    error("Move \"$(move)\" not recognized.")
end


"""
Randomly generates a move proposal from the function `sample_move_proposal`.
Sometimes the move sampled in `sample_move_proposal` will not have any valid
parameters, so this function repeatedly samples until a move with valid
parameters is selected. Also, makes sure the proposed move is not present in
`moves_attempted`, and appends the returned move to this vector.
"""
function generate_move_proposal(N::HybridNetwork, moves_attempted::Vector, hmax::Int, rng::TaskLocalRNG)
    retries::Int = 0
    move, params = sample_move_proposal(N, hmax, rng)

    while params === nothing || already_attempted(moves_attempted, move, params)
        move, params = sample_move_proposal(N, hmax, rng)
        retries += 1
        if retries >= 1e6
            error("Could not find any valid move proposals after 1e6 attempts.")
        end
    end

    return (move, params)
end


"""
Helper function that determines whether the move `move` with parameters `params` has
already been attempted (i.e. is stored in the vector `moves_attempted`).
"""
function already_attempted(moves_attempted::Vector, move::Symbol, params::Tuple)::Bool
    for (amove, aparams) in moves_attempted
        move == amove || continue
        all_params_match::Bool = true
        for (aparam, param) in zip(aparams, params)
            if !(typeof(aparam) <: typeof(param))
                all_params_match = false
                break
            end
            if aparam.number != param.number || (typeof(aparam) <: Node && aparam.name != param.name)
                all_params_match = false
                break
            end
        end
        all_params_match && return true
        break
    end
    return false
end


"""
Randomly samples a move to generate a new topology from `N`.
"""
function sample_move_proposal(N::HybridNetwork, hmax::Int, rng::TaskLocalRNG)

    if N.numhybrids < hmax && rand(rng) < 0.05
        @debug "SELECTED: add_random_hybrid!"
        return (:add_hybrid, sample_add_hybrid_parameters(N, rng))
    end

    # If net has 0 hybrids, we can only do rNNI(1) or rSPR moves
    if N.numhybrids == 0
        # PROBABILITY OF EACH MOVE:
        # rNNI(1):  70%
        # rSPR:     30%

        r = rand(rng)
        if r < 0.7
            return (:rNNI1, sample_rNNI_parameters(N, 1, rng))
        else
            return (:rSPR, sample_rSPR_parameters(N, rng))
        end
    end


    # PROBABILITY OF EACH MOVE:
    # rNNI(1):      15%
    # rNNI(2):      0%  NEVER accepted!
    # rNNI(3):      0%  NEVER accepted!
    # rNNI(4):      10%
    # rSPR:         5%
    # origin:       10%
    # target:       10%
    # loc origin:   20%
    # loc target:   15%
    # fliphybrid:   15%
    probs = [0.15, 0.0, 0.0, 0.1, 0.05, 0.1, 0.1, 0.2, 0.15, 0.15]

    r = rand(rng)
    if r < sum(probs[1:1])
        return (:rNNI1, sample_rNNI_parameters(N, 1, rng))
    elseif r < sum(probs[1:2])
        return (:rNNI2, sample_rNNI_parameters(N, 2, rng))
    elseif r < sum(probs[1:3])
        return (:rNNI3, sample_rNNI_parameters(N, 3, rng))
    elseif r < sum(probs[1:4])
        return (:rNNI4, sample_rNNI_parameters(N, 4, rng))
    elseif r < sum(probs[1:5])
        return (:rSPR, sample_rSPR_parameters(N, rng))
    elseif r < sum(probs[1:6])
        return (:retic_origin, sample_move_reticulate_origin_parameters(N, rng))
    elseif r < sum(probs[1:7])
        return (:retic_target, sample_move_reticulate_target_parameters(N, rng))
    elseif r < sum(probs[1:8])
        return (:retic_origin_local, sample_move_reticulate_origin_local_parameters(N, rng))
    elseif r < sum(probs[1:9])
        return (:retic_target_local, sample_move_reticulate_target_local_parameters(N, rng))
    else
        return (:flip_hybrid, sample_flip_hybrid_parameters(N, rng))
    end

end


