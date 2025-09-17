using PhyloNetworks, SNaQ, DataFrames, Random, Distributed
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
true_logPL = compute_loss(net, q)   # -0.00696866191062832
#############



#### FOR LOOP HERE
for_loop_rt = @elapsed for_loop_best, for_loop_nets, _ = multi_search(deepcopy(gts[1]), q, 1; seed=42, runs=10, maxeval=10000, maxequivPLs=500)
[hardwiredclusterdistance(for_loop_nets[j], net, false) for j in eachindex(for_loop_nets)]

@everywhere include("../../src/gradient_optimization/search_API.jl")
@everywhere __precompile__()
pmap_rt = @elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; seed=0, runs=10, maxeval=10000, maxequivPLs=500)
println([hardwiredclusterdistance(pmap_nets[j], net, false) for j in eachindex(pmap_nets)])


@everywhere include("../../src/gradient_optimization/search_API.jl")

@elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; propQuartets=0.5, seed=0, runs=5, maxeval=10000, maxequivPLs=500)
@elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; propQuartets=1.0, seed=0, runs=5, maxeval=10000, maxequivPLs=500)


println([hardwiredclusterdistance(pmap_nets[j], net, false) for j in eachindex(pmap_nets)])





#[0, 0, 4, 4, 6, 6, 6, 8, 8, 8]

# 3 proc, 4 thread:  31.31s
# 3 proc, 1 thread:  73.44s
# 1 proc, 4 thread:  69.51s
# 1 proc, 1 thread: 157.14s




# Correct h
Random.seed!(42)
opt_rt = @elapsed best_net, est_nets, logPLs = multi_search(deepcopy(gts[1]), q, 1; runs=10, maxeval=10000, maxequivPLs=500)
@info [hardwiredclusterdistance(n, net, false) for n in est_nets]

snaq_rt = @elapsed snaq_net = snaq!(gts[1], df, hmax=1, runs=10)
@info hardwiredclusterdistance(snaq_net, net, false)
snaq_rt, opt_rt

all_snaq_nets = []
for line in readlines("snaq.out")[6:(length(readlines("snaq.out"))-2)]
    line = line[2:length(line)]
    line = split(line, ", with -loglik")[1]
    push!(all_snaq_nets, readnewick(line))
end
[hardwiredclusterdistance(snet, net, false) for snet in all_snaq_nets]

@info 2*snaq_rt / opt_rt










