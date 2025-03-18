using PhyloNetworks, Plots, SNaQ, DataFrames, Random
include("../../src/network_moves/rNNI_moves.jl")
include("../../src/network_moves/rSPR_moves.jl")
include("../../src/gradient_optimization/opt_API.jl")
include("../../src/gradient_optimization/search_API.jl")
ENV["JULIA_DEBUG"] = Main
ENV["JULIA_DEBUG"] = nothing

using PhyloCoalSimulations
using PhyloPlots

### SETUP ###
net, gts = readnewick(joinpath(@__DIR__, "n1.netfile")), readmultinewick(joinpath(@__DIR__, "n1_10kgts.treefile")); 
net.isrooted = false;
q, t = countquartetsintrees(gts, showprogressbar=false);
df = readtableCF(DataFrame(tablequartetCF(q, t)));
tre0 = readnewick(writenewick(gts[1]));
#perform_random_rNNI!(tre0);
true_logPL = compute_loss(net, q)
#############


# Correct h
Random.seed!(42)
opt_rt = @elapsed best_net, est_nets, logPLs = multi_search(deepcopy(gts[1]), q, 1; runs=10, maxeval=10000, maxequivPLs=500)
@info [hardwiredClusterDistance(n, net, false) for n in est_nets]

snaq_rt = @elapsed snaq_net = snaq!(gts[1], df, hmax=1, runs=10)
@info hardwiredClusterDistance(snaq_net, net, false)
snaq_rt

all_snaq_nets = []
for line in readlines("snaq.out")[6:(length(readlines("snaq.out"))-2)]
    line = line[2:length(line)]
    line = split(line, ", with -loglik")[1]
    push!(all_snaq_nets, readnewick(line))
end
[hardwiredClusterDistance(snet, net, false) for snet in all_snaq_nets]

@info 2*snaq_rt / opt_rt







times = [];
accs = [];
max_equivs = [];
max_eval = [];
for maxequiv in [10, 50] # [10, 50, 100, 250]
    for maxeval in [1, 3, 5, 10, 25, 50]
        for _ = 1:10
            push!(max_equivs, maxequiv)
            opt_rt = @elapsed best_net, _ = multi_search(deepcopy(gts[1]), q, 1; runs=1, maxeval=10000, maxequivPLs=maxequiv)
            push!(times, opt_rt)
            push!(accs, hardwiredClusterDistance(best_net, net, false))
        end
    end
end

using Plots




