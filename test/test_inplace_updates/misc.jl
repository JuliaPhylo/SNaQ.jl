using PhyloNetworks, PhyloCoalSimulations, Random, BenchmarkTools
include("../../src/network_moves/add_remove_retic.jl")


function generate_tree(n::Int, seed::Int=42)
    n > 0 || error("n must be >0 (n=$(n))")
    Random.seed!(seed)

    newick = "(" * join(["t$(j)" for j=1:n], ",") * ");"
    tre0 = readnewick(newick)
    for E in tre0.edge
        E.length = 0.0
    end
    return simulatecoalescent(tre0, 1, 1)[1]
end


function generate_net(n::Int, h::Int, seed::Int=42)
    (n > 0 && h > 0) || error("n and h must be >0 (n=$(n), h=$(h))")
    rng = Random.seed!(seed)

    net = generate_tree(n, seed)
    for j = 1:h
        add_random_hybrid!(net, rng)
    end
    return net
end

function opt_benchmark(ntaxa::Int=15, nhyb::Int=3, seed::Int=42)
    Random.seed!(seed)
    net = generate_net(ntaxa, nhyb, seed)
    gts = simulatecoalescent(net, 1000, 1);
    q, t = countquartetsintrees(gts, showprogressbar=false);
    blocks, _ = find_quartet_equations(net);

    return @benchmark optimize_bls!(net, blocks, q)
end

function opt_profile(ntaxa::Int=15, nhyb::Int=3, seed::Int=42; nrep::Int=1)
    Random.seed!(seed)
    net = generate_net(ntaxa, nhyb, seed)
    gts = simulatecoalescent(net, 1000, 1);
    q, t = countquartetsintrees(gts, showprogressbar=false);
    blocks, _ = find_quartet_equations(net);

    start_time = time()
    @profview for _ = 1:nrep optimize_bls!(net, blocks, q) end
    @info "Profiling took $(round(time() - start_time, digits=4)) seconds."
end


function eqn_profile(ntaxa::Int=15, nhyb::Int=3, seed::Int=42; nrep::Int=1)
    Random.seed!(42)
    net = generate_net(ntaxa, nhyb, seed)

    start_time = time()
    @profview for _ = 1:nrep find_quartet_equations(net) end
    @info "Profiling took $(round(time() - start_time, digits=4)) seconds."
end

