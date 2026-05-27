"""
    multisearch(
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
    of topologies. If a vector, must be length eor of length `runs`.
- `q::Union{DataCF, AbstractArray{Float64}}`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `runs::Int=10`: Number of independent search runs.
- `seed::Int=42`: Random seed for reproducibility.
- `restrictions::Fuction=defaultrestrictions()`: Function that takes a `HybridNetwork` as
    its only argument and returns a `Bool`. Only networks that return `true` to this function
    will be considered during search.
- `kwargs...`: Additional keyword arguments passed to the [`search`](@ref) function.

# Returns
- `best_network::HybridNetwork`: The network with the best (lowest) negative log pseudo-likelihood.
- `all_networks::Vector{HybridNetwork}`: All networks from the runs, sorted by score.
- `all_scores::Vector{Float64}`: Negative log pseudo-likelihood scores for each network.

# Notes
- This function uses distributed computing to perform searches in parallel.
- It returns the best network found across all runs.
"""
function multisearch(
    N::Union{HybridNetwork, AbstractVector{HybridNetwork}},
    q::Union{DataCF, AbstractMatrix{Float64}},
    hmax::Int;
    # Basic arguments
    runs::Int=10,
    seed::Int=rand(Int),
    logprefix::String="",
    filename::String="snaq",
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
        q = gatherCFmatrix(q)
    end

    # Generate per-run seeds
    Random.seed!(seed)
    run_seeds = abs.(rand(Int, runs) .% 100000)

    # Log run details
    if filename != ""
        open("$(filename).log", "w+") do f end
        restrictionmsg = restrictions == defaultrestrictions() ? "default restrictions" :
            restrictions == restrictgallednetwork() ? "galled networks" :
            restrictions == restrictgalledtree() ? "galled trees" :
            restrictions == restrictrootedtreechild() ? "rooted tree child" :
            restrictions == restrictweaklytreechild() ? "weakly tree child" :
            restrictions == restrictstronglytreechild() ? "strongly tree child" : "custom restrictions"
        
        logmessage(filename, """
        Beginning network optimization using SNaQ.jl across $runs runs with the following parameters:
            hmax = $hmax,
            seed = $seed,
            outgroup = $outgroup,
            restrictions = $restrictionmsg
        Root name for log files: $filename (absolute path $(abspath(filename)))
        Currently utilizing $(nprocs()) processor$(nprocs() > 1 ? "s" : "") and $(Threads.nthreads()) thread$(Threads.nthreads() > 1 ? "s" : "").
        """)
    end

    # Make sure the log files can be created
    runs_path::String = ""
    if filename != ""
        runs_path = string(filename, "_runs/")
        mkpath(runs_path)
    end
    filenames = [filename == "" ? "" : "$(runs_path)run$(j)" for j = 1:runs]

    # Do the runs distributed
    starttime = time()
    all_nets = pmap(
        j -> search(
            length(Ns) == 1 ? Ns[1] : Ns[j],
            q, hmax; seed = run_seeds[j], restrictions=restrictions,
            logfile = logprefix == "" ? "" : "$(logprefix)$(j)",
            filename = filenames[j],
            outgroup=outgroup, kwargs...
        ),
        1:runs
    )
    elapsed = timeelapsed(time() - starttime)

    # Consolidate return data
    sort_idx = sortperm(loglik.(all_nets), rev=true)
    bestnet = all_nets[sort_idx[1]]

    # Log results
    logmessage(filename, """
    Finished optimizing topology at $(currenttime()) after $(elapsed).
    Optimal network: $(writenewick(bestnet, round=true))
    Optimal -loglik: $(-loglik(bestnet))
    To view all $runs inferred networks and their associated -loglik scores, see $(filename).out ($(abspath("$(filename).out")))""")

    open("$(filename).out", "w+") do f
        print(f,
            """
            $(writenewick(bestnet)) -Ploglik = $(-loglik(bestnet))
             Elapsed time: $(elapsed), $(runs) attempted runs
            
            -----------------------------------
            List of estimated networks for all runs (sorted by log-pseudolik; the smaller, the better):
            """
        )
        for j in sort_idx
            println(f, " $(writenewick(all_nets[j])), with -loglik $(loglik(all_nets[j]))")
        end
        println(f, "-----------------------------------")
    end

    open("$(filename).networks", "w+") do f
        for i in sort_idx
            write(f, "$(writenewick(all_nets[i])), with -loglik $(loglik(all_nets[i]))")
            if i == 1
                write(f, " (best network found, remaining sorted by log-pseudolik; the smaller, the better)")
            end
            write(f, "\n")
        end
    end

    # Clean up: the edges above roots have leftover values in them right now -
    #           we can't actually infer the lengths of those edges, so we clean
    #           those up here.
    for n in all_nets
        for L in n.leaf
            getparentedge(L).length = -1
        end
    end

    # Return
    return bestnet, all_nets[sort_idx]
end


"""
Verifies that starting topology(ies) `N` are ready to be optimized.
Modifies the network(s) `N` in-place.
"""
function verifystartingtopologies!(N::Union{HybridNetwork, AbstractVector{HybridNetwork}}, outgroup::String, restrictions::Function)::Vector{HybridNetwork}
    # Copy the input networks
    Ns::Vector{HybridNetwork} = typeof(N) <: HybridNetwork ? [deepcopynetwork(N)] : [deepcopynetwork(n) for n in N]
    for (j, n) in enumerate(Ns)
        # Split multifurcations
        if any(n -> length(n.edge) > 3, Ns[j].node)
            @warn "Input network #$(j) has a polytomy. SNaQ only infers binary networks, so this will be automatically resolved before inference."
            
            iters = 0
            while true
                iters += 1
                if iters > N.numnodes
                    error("Got stuck in an infinite loop while resolving polytomies. Please report this bug with your input tree(s) on GitHub.")
                end
                multi = findfirst(n -> length(n.edge) > 3, Ns[j].node)
                if isnothing(multi) break end
                PhyloNetworks.resolvetreepolytomy!(Ns[j], Ns[j].node[multi])
            end
        end

        # Prep data
        semidirectnetwork!(Ns[j]);

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
function logtext(logfile::String, msg::String)
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
function logmoves(logfile::String, moves_prop::Dict, moves_acc::Dict, moves_PL::Dict)
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

    logtext(logfile, msg)
end

logmessage(filename::String, msg::String) = remotecall_fetch(writelogmessage, 1, filename, msg)
currenttime() = Dates.format(now(), "HH:MM:SS yyyy-mm-dd")

function writelogmessage(filename::String, msg::String)
    open("$(filename).log", "a+") do f
        println(f, msg)
    end
end

function timeelapsed(totaltime::Float64)::String
    seconds::Int = Int(round(totaltime))
    minutes::Int = (seconds ÷ 60) % 60
    hours::Int   = (seconds ÷ 3600) % 24
    days::Int    = seconds ÷ 86400
    seconds      = seconds % 60
    if days > 0
        return "$days days, $hours hours, $minutes minutes and $seconds seconds"
    elseif hours > 0
        return "$hours hours, $minutes minutes and $seconds seconds"
    elseif minutes > 0
        return "$minutes minutes and $seconds seconds"
    else
        return "$seconds seconds"
    end
end



"""
Performs a single search for the optimal network topology with gradient-based optimization
of branch lengths and inheritance probabilities.

# Arguments
- `N::HybridNetwork`: The starting network topology.
- `q`: Observed quartet concordance factors.
- `hmax::Int`: Maximum number of hybridization events allowed.

# Optional Arguments
- `restrictions::Function=defaultrestrictions()`: Function to enforce restrictions on the proposed networks.
- `qinfTest::Bool`: whether to test for uninformative quartets (CFs near [1/3, 1/3, 1/3]) and
    exclude those quartets when performing the network search.
- `qtolAbs::Float64=1e-4`: the absolute tolerance used to detect uninformative quartets.
- `propQR::Float64=0.0`: probability at a given search iteration of utilizing weighted random
    sampling to (i) sample a poorly fitting quartet, (ii) sample an edge spanned by that
    quartet in the network, and finally (iii) only sample moves that include this edge.
- `α::Real=Inf`: Dirichlet parameter for gene tree heterogeneity model.
- `propQuartets::Real=1.0`: Proportion of quartets to use during optimization.
- `preopt::Bool=false`: Whether to perform a pre-optimization step.
- `probST::Real=0.3`: Probability of performing a subtree move before searching.
- `maxeval::Int=Int(1e8)`: Maximum number of evaluations.
- `maxequivPLs::Int=1500`: Maximum number of equivalent pseudo-likelihood scores to consider.
- `opt_maxeval::Int=30`: Maximum evaluations for optimization.
- `seed::Int=abs(rand(Int) % 100000)`: Random seed for reproducibility.
- `verbose::Bool=false`: Whether to print verbose output.
- `logfile::String=""`: File to log detailed progress (used for debugging, but can also be used to examine convergence).

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
    preopt::Bool=true,
    probST::Real=0.3,
    probQR::Float64=0.0,
    maxeval::Int=Int(1e8),
    maxequivPLs::Int=1500,
    liktolAbs::Float64=1e-8,
    liktolRel::Float64=1e-4,
    opt_maxeval::Int=30,
    seed::Int=abs(rand(Int) % 100000),
    verbose::Bool=false,
    logfile::String="",
    filename::String="",
    outgroup::String="none",
    qinfTest::Bool=false,
    qtolAbs::Float64=1e-4,
    optargs...
)
    # Parameter enforcement
    maxeval > 0 || error("maxeval must be > 0 (maxeval = $(maxeval)).")
    maxequivPLs > 0 || error("maxequivPLs must be > 0 (maxequivPLs = $(maxequivPLs)).")
    0 ≤ α ≤ Inf || error("α must be in range [1, ∞] (α = $(α))")
    0 < propQuartets ≤ 1 || error("propQuartets must be in range (0, 1] (propQuartets = $(propQuartets))")
    0 ≤ probQR ≤ 1 || error("probQR must be in range [0, 1] (probQR = $(probQR))")
    0 ≤ probST ≤ 1 || error("probST must be in range [0, 1] (probST = $(probST))")
    outgroup == "none" || any(l -> l.name == outgroup, N.leaf) || error("No taxa in N have taxa name $(outgroup) (outgroup name)")
    qtolAbs ≥ 0.0 || error("qtolAbs must be ≥ 0.0 (qtolAbs = $qtolAbs)")
    0 ≤ probQR ≤ 1 || error("probQR must be in range [0, 1] (probQR = $probQR)")

    # Initial logging message
    starttime = time()
    open(string(filename, ".log"), "w+") do f end # clear any pre-existing text in the log file
    logmessage(filename, """
    BEGIN: search with seed $(seed) at $(currenttime())
           starting topology: $(writenewick(N, round=true))""")

    # Convert q to a Matrix if it is a DataCF
    if typeof(q) <: DataCF
        q = gatherCFmatrix(q)
    end

    # Set the seed
    rng = Random.seed!(seed)

    # N = readnewick(writenewick(N));
    N = deepcopynetwork(N);
    for node in N.node
        if !node.leaf && !node.hybrid
            node.name = ""
        end
    end
    semidirectnetwork!(N)
    restrictions(N) || error("N does not meet restrictions IMMEDIATELY")

    if rand(rng) < probST
        found_different_net::Bool = false
        for j = 1:10
            try
                performrNNI1!(N, samplerNNIparameters(N, 1, rng)...);
                if restrictions(N)
                    found_different_net = true
                    break
                end
            catch
            finally
                if !found_different_net
                    # N = readnewick(writenewick(N));
                    N = deepcopynetwork(N);
                    semidirectnetwork!(N)
                end
            end
        end
        if !found_different_net
            @warn "Initial probST move led to a network that did not meet the given restrictions. Using the provided network instead."
        end
    end

    # Data used throughout the optimization process
    local q_idxs::Vector{Int64}
    if qinfTest
        informative = trues(size(q, 1))
        if qinfTest
            for (i, row) in enumerate(eachrow(q))
                informative[i] = isquartetinformative(row, qtolAbs)
            end
        end
        q_idxs = sampleqindices(N, propQuartets, informative, rng)
    elseif propQuartets == 1.0
        q_idxs = collect(1:nchoose4taxalength(N))
    else
        q_idxs = sampleqindices(N, propQuartets, rng)
    end
    logPLs::Array{Float64} = Array{Float64}(undef, maxeval)
    neq = findquartetequations(N, q_idxs);
    N_eqns::Vector{QuartetData} = neq[1];
    CFΔs::Vector{Float64} = []  # used when probQR != 0.0, computed WHEN NEEDED, so init'd to []
    unchanged_iters = 0

    # Pre-optimizing the network's parameters
    if preopt
        @debug "Pre-optimizing"
        optimize!(N, N_eqns, q[q_idxs, :], α; maxeval=max(opt_maxeval, 500), optargs...)
        restrictions(N) || error("N does not meet restrictions after preopt")
        logPLs[1] = loglik(N)
    else
        logPLs[1] = computeloss(N_eqns, gatherparams(N), q[q_idxs, :], α)
    end

    moves_attempted = [];   # Vector of Tuples: (<move name>, <move parameters (i.e. nodes/edges)>)
    moves_proposed = Dict{Symbol,Int}()
    moves_accepted = Dict{Symbol,Int}()
    moves_logPL = Dict{Symbol,Vector{Float64}}()
    last_move = :nothing

    logtext(logfile, "Entering main loop with -logPL = $(logPLs[1])")
    for j = 2:maxeval
        if j % 100 == 0
            logmoves(logfile, moves_proposed, moves_accepted, moves_logPL)
        end

        verbose && print("\rIteration $(j)/$(maxeval) - in a row=$(unchanged_iters)/$(maxequivPLs)              ")

        # 1. Propose a new topology
        @debug "Current: $(writenewick(N, round=true))"
        #Nprime = readnewick(writenewick(N));
        Nprime = deepcopynetwork(N);

        prop_move, prop_params = generatemoveproposal(Nprime, N_eqns, moves_attempted, hmax, probQR, q[q_idxs, :], CFΔs, rng, α)
        last_move = prop_move
        applymove!(Nprime, prop_move, prop_params)
        @debug "Proposed move: $(prop_move), parameters: $(prop_params)"
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
                    logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (cannot reroot at outgroup)")
                    logPLs[j] = logPLs[j-1]
                    continue
                else
                    rethrow(e)
                end
            end
        end

        # 2.3 After removing some edges above, the root may have 2 edge now instead of 3 - we fix that here
        semidirectnetwork!(Nprime)

        # 3. Immediately throw away networks that don't meet restrictions 
        if !restrictions(Nprime)
            @debug "Nprime does not meet restrictions - skipping."
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) (restrictions not met)")
            logPLs[j] = logPLs[j-1]
            continue
        end

        # Check whether we can do in-place updates here.
        # We CANNOT do inplace updates if:
        # 1. the number of hybrids changes, OR
        # 2. the number of optimization parameters in the network changed
        Nprime_np::Int = gatheroptimizationinfo(Nprime)[1]
        N_np::Int = gatheroptimizationinfo(N)[1]
        cannot_do_inplace::Bool = N.numhybrids != Nprime.numhybrids || N_np != Nprime_np

        # 4. Optimize branch lengths and compute logPL
        Nprime_logPL, Nprime_eqns = optimizetopology!(
            Nprime, N_eqns, prop_move, prop_params, q, q_idxs,
            opt_maxeval, cannot_do_inplace, rng, α; optargs...
        )
        Nprime_logPL == -Inf && error("Nprime_logPL is -Inf?? newick: $(writenewick(Nprime, round=true))\nold network: $(writenewick(N, round=true))\nprop move: $(prop_move)\nprop params: $(prop_params)")
        # computeloss(Nprime, q) == Nprime_logPL || error("LOGPLS NOT EQUAL AFTER MOVE $(prop_move)")

        # 5. Accept / reject
        isnan(Nprime_logPL) && error("""
            Nprime_logPL = $(Nprime_logPL)
            $(writenewick(Nprime, round=true))
        """)
        if Nprime_logPL - logPLs[j-1] > liktolAbs && (logPLs[j-1] - Nprime_logPL) / logPLs[j-1] > liktolRel
            # Update current topology info
            N = Nprime
            N_eqns = Nprime_eqns
            logPLs[j] = Nprime_logPL
            CFΔs = []
            moves_accepted[prop_move] += 1
            push!(moves_logPL[prop_move], logPLs[j] - logPLs[j-1])

            # Log acceptance
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) ACCEPTED $(prop_move), new -logPL=$(round(logPLs[j], digits=6))")

            # Update tracking vars
            unchanged_iters = 0
            moves_attempted = []
        else
            logPLs[j] = logPLs[j-1]
            unchanged_iters += 1

            # Log rejection and reason
            logtext(logfile, "Iteration $(j) (N.h=$(N.numhybrids)), in a row = $(unchanged_iters)/$(maxequivPLs) REJECTED $(prop_move) ($(round(Nprime_logPL, digits=3)) < $(round(logPLs[j], digits=3)))")
        end

        # Early stopping checks
        if unchanged_iters > maxequivPLs
            @debug "stopping early after $(j) iterations"
            logPLs = logPLs[1:j]
            break
        end
    end
    loglik!(N, logPLs[length(logPLs)])

    if propQuartets != 1.0
        logmessage(filename, "Re-optimizing branch lengths with ALL quartets.")
        if typeof(q) <: DataCF
            loglik!(N, optimize!(N, gatherCFmatrix(q)))
        else
            loglik!(N, optimize!(N, q))
        end
        logmessage(filename, "END propQuartets<1.0 post-search parameter optimization: found minimizer topology with -loglik=$(round(loglik(N), digits=5))")
    end

    # Remove internal node names that are not hybrids
    for node in N.node
        if !node.leaf && !node.hybrid
            node.name = ""
        end
    end

    # Rename hybrids to be from 1-H
    for (iH, H) in enumerate(N.hybrid)
        H.name = "H$iH"
    end

    logtext(logfile, "Search complete at $(currenttime()).\n\n")
    logmoves(logfile, moves_proposed, moves_accepted, moves_logPL)

    logmessage(filename, "END: search with seed $(seed) after $(timeelapsed(time() - starttime)). -Ploglik=$(-loglik(N))")
    logmessage(filename, writenewick(N))
    return N
end


"""
Checks whether a given quartet is informative based off of observed CFs.
Informative here is defined as any two entries in the quartet's observed
CFs having absolute difference greater than `atol`.
"""
function isquartetinformative(ocfrow::AbstractVector{Float64}, atol::Float64)
    for i = 1:2
        for j = (i+1):3
            if abs(ocfrow[i] - ocfrow[j]) > atol
                return true
            end
        end
    end
    return false
end


"""
Applies the move `move` on parameters `params` to network `N`.
"""
function applymove!(N::HybridNetwork, move::Symbol, params::Tuple)
    if move == :addhybrid
        return addhybrid!(N, params...)
    elseif move == :rNNI1
        return performrNNI1!(N, params...)
    elseif move == :rNNI2
        return performrNNI2!(N, params...)
    elseif move == :rNNI3
        return performrNNI3!(N, params...)
    elseif move == :rNNI4
        return performrNNI4!(N, params...)
    elseif move == :retic_origin || move == :retic_origin_local
        return movereticulateorigin!(N, params...)
    elseif move == :retic_target || move == :retic_target_local
        return movereticulatetarget!(N, params...)
    elseif move == :rSPR
        return performrSPR!(N, params...)
    elseif move == :flip_hybrid
        return fliphybrid!(N, params[1])
    end

    error("Move \"$(move)\" not recognized.")
end


"""
Randomly generates a move proposal from the function `samplemoveproposal`.
Sometimes the move sampled in `samplemoveproposal` will not have any valid
parameters, so this function repeatedly samples until a move with valid
parameters is selected. Also, makes sure the proposed move is not present in
`moves_attempted`, and appends the returned move to this vector.
"""
function generatemoveproposal(Nprime::HybridNetwork, N_eqns::Vector{QuartetData}, moves_attempted::Vector, hmax::Int, probQR::Float64, Q::Matrix{Float64}, CFΔs::Vector{Float64}, rng::TaskLocalRNG, α::Float64)::Tuple{Symbol,Any}
    required_edge::Union{Edge, Nothing} = nothing
    if probQR > 0.0 && rand(rng) <= probQR
        required_edge = sampleprobQRedge(Nrpime, N_eqns, Q, CFΔs, rng, α)
    end

    validmove(mv::Symbol, pars) = !isnothing(pars) &&
        !alreadyattempted(moves_attempted, mv, pars) &&
        (isnothing(required_edge) || any(p -> p == required_edge, pars))

    retries::Int = 0
    move, params = samplemoveproposal(Nprime, hmax, rng)

    while !validmove(move, params)
        move, params = samplemoveproposal(Nprime, hmax, rng)
        retries += 1
        if retries >= 1e6
            error("Could not find any valid move proposals after 1e6 attempts.")
        end
    end

    return (move, params)
end
generatemoveproposal(Nprime::HybridNetwork, ma::Vector, hmax::Int, rng::TaskLocalRNG) =
    generatemoveproposal(Nprime, Vector{QuartetData}([]), ma, hmax, 0.0, zeros(0, 0), zeros(0), rng, Inf)


"""
Samples an edge from the network `N` with weights stored in `CFΔs`. If `CFΔs` is empty
(it starts empty and is reset to [] whenever the search finds a better network), these
weights are computed.
"""
function sampleprobQRedge(N::HybridNetwork, eqns::Vector{QuartetData}, Q::Matrix{Float64}, CFΔs::Vector{Float64}, rng::TaskLocalRNG, α::Float64)::Edge
    if length(CFΔs) == 0
        params = gatherparams(N);
        for (eqn, (ocf1, ocf2, ocf3)) in zip(eqns, eachrow(Q))
            ecf1, ecf2, ecf3 = computeexpectedCF(eqn, params, α)
            push!(CFΔs, abs(ecf1 - ocf1) + abs(ecf2 - ocf2) + abs(ecf3 - ocf3))
        end
    end

    # Queue approach to run through every contributing equation in the sampled
    # CF's equation to make sure we select ALL edges that relate to this quartet
    quartet = sample(rng, 1:length(eqns), Weights(CFΔs))
    Q = [quartet.eqn]
    edges = []
    while length(Q) > 0
        curr = Q[1]
        deleteat!(Q, 1)

        append!(edges, curr.coal_edges)
        append!(Q, curr.divisions)
    end

    return sample(rng, unique(edges))
end


"""
Helper function that determines whether the move `move` with parameters `params` has
already been attempted (i.e. is stored in the vector `moves_attempted`).
"""
function alreadyattempted(moves_attempted::Vector, move::Symbol, params::Tuple)::Bool
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
function samplemoveproposal(N::HybridNetwork, hmax::Int, rng::TaskLocalRNG)::Tuple{Symbol,Any}
    if N.numhybrids < hmax && rand(rng) < 0.05
        @debug "SELECTED: add_random_hybrid!"
        return (:addhybrid, sampleaddhybridparameters(N, rng))
    end

    # If net has 0 hybrids, we can only do rNNI(1) or rSPR moves
    if N.numhybrids == 0
        # PROBABILITY OF EACH MOVE:
        # rNNI(1):  70%
        # rSPR:     30%

        r = rand(rng)
        if r < 0.7
            return (:rNNI1, samplerNNIparameters(N, 1, rng))
        else
            return (:rSPR, samplerSPRparameters(N, rng))
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
    # local origin: 20%
    # local target: 15%
    # fliphybrid:   15%
    probs = [0.15, 0.0, 0.0, 0.1, 0.05, 0.1, 0.1, 0.2, 0.15, 0.15]
    cumprobs = cumsum(probs)

    r = rand(rng)
    if r <= sum(cumprobs[1])
        return (:rNNI1, samplerNNIparameters(N, 1, rng))
    elseif r <= sum(cumprobs[2])
        return (:rNNI2, samplerNNIparameters(N, 2, rng))
    elseif r <= sum(cumprobs[3])
        return (:rNNI3, samplerNNIparameters(N, 3, rng))
    elseif r <= sum(cumprobs[4])
        return (:rNNI4, samplerNNIparameters(N, 4, rng))
    elseif r <= sum(cumprobs[5])
        return (:rSPR, samplerSPRparameters(N, rng))
    elseif r <= sum(cumprobs[6])
        return (:retic_origin, samplemovereticulateoriginparameters(N, rng))
    elseif r <= sum(cumprobs[7])
        return (:retic_target, samplemovereticulatetargetparameters(N, rng))
    elseif r <= sum(cumprobs[8])
        return (:retic_origin_local, samplemovereticulateoriginlocalparameters(N, rng))
    elseif r <= sum(cumprobs[9])
        return (:retic_target_local, samplemovereticulatetargetlocalparameters(N, rng))
    else
        return (:flip_hybrid, samplefliphybridparameters(N, rng))
    end

end


