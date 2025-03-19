using PhyloNetworks, Plots, SNaQ, DataFrames, Random, Distributed
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
[hardwiredClusterDistance(for_loop_nets[j], net, false) for j in eachindex(for_loop_nets)]

@everywhere include("../../src/gradient_optimization/search_API.jl")
@everywhere __precompile__()
pmap_rt = @elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; seed=0, runs=10, maxeval=10000, maxequivPLs=500)
println([hardwiredClusterDistance(pmap_nets[j], net, false) for j in eachindex(pmap_nets)])


@everywhere include("../../src/gradient_optimization/search_API.jl")

@elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; propQuartets=0.5, seed=0, runs=5, maxeval=10000, maxequivPLs=500)
@elapsed pmap_best, pmap_nets, _ = multi_search(deepcopy(gts[1]), q, 1; propQuartets=1.0, seed=0, runs=5, maxeval=10000, maxequivPLs=500)


println([hardwiredClusterDistance(pmap_nets[j], net, false) for j in eachindex(pmap_nets)])





#[0, 0, 4, 4, 6, 6, 6, 8, 8, 8]

# 3 proc, 4 thread:  31.31s
# 3 proc, 1 thread:  73.44s
# 1 proc, 4 thread:  69.51s
# 1 proc, 1 thread: 157.14s




# Correct h
Random.seed!(42)
opt_rt = @elapsed best_net, est_nets, logPLs = multi_search(deepcopy(gts[1]), q, 1; runs=10, maxeval=10000, maxequivPLs=500)
@info [hardwiredClusterDistance(n, net, false) for n in est_nets]

snaq_rt = @elapsed snaq_net = snaq!(gts[1], df, hmax=1, runs=10)
@info hardwiredClusterDistance(snaq_net, net, false)
snaq_rt, opt_rt

all_snaq_nets = []
for line in readlines("snaq.out")[6:(length(readlines("snaq.out"))-2)]
    line = line[2:length(line)]
    line = split(line, ", with -loglik")[1]
    push!(all_snaq_nets, readnewick(line))
end
[hardwiredClusterDistance(snet, net, false) for snet in all_snaq_nets]

@info 2*snaq_rt / opt_rt







using DataFrames, CSV
df = CSV.read(joinpath(@__DIR__, "10run_param_perf.csv"), DataFrame)

# This if CSV is empty
# df = DataFrame(
#     times=Float64[], accs=Int[], mean_accs=Float64[], median_accs=Float64[],
#     max_equivs=Int[], max_evals=Int[], ntaxa=Int[], nhyb=Int[]
# )


include("../test_inplace_updates/misc.jl")
nruns = 10
while true
    for ntaxa in [8, 10, 12, 15]
        for nhyb in [1, 2, 3]
            net = generate_net(ntaxa, nhyb, abs(rand(Int)))
            rootonedge!(net, getroot(net).edge[1])
            gts = simulatecoalescent(net, 1000, 1);
            q, t = countquartetsintrees(gts; showprogressbar=false)
            semidirect_network!(net)
            println("n$(ntaxa)h$(nhyb)")

            for maxequiv in [100, 500, 1000]
                for maxeval in [10, 25, 50]
                    print("\rmaxPL=$(maxequiv), maxeval=$(maxeval)        ")

                    opt_rt = @elapsed best_net, all_nets, _ = multi_search(deepcopy(gts[1]), q, 1; runs=nruns, maxeval=10000, maxequivPLs=maxequiv, opt_maxeval=maxeval)
                    hwcds = [hardwiredClusterDistance(opt_net, net, false) for opt_net in all_nets]
                    
                    push!(df, [
                        opt_rt, hardwiredClusterDistance(best_net, net, false),
                        mean(hwcds), median(hwcds), maxequiv,
                        maxeval, ntaxa, nhyb
                    ])
                    CSV.write(joinpath(@__DIR__, "10run_param_perf.csv"), df)

                    #p = Gadfly.plot(df, x=:max_equivs, y=:times, color=:max_evals, Geom.point)
                    # stack_df = stack(df, [:times, :accs], variable_name="metric", value_name="value")
                    # p = Gadfly.plot(stack_df, x=:max_evals, y=:value, xgroup=:max_equivs, ygroup=:metric, Geom.subplot_grid(Geom.boxplot))
                    # display(p)
                end
            end
            println("")
        end
    end
end





