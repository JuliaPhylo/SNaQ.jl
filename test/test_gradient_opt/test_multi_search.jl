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







using DataFrames, Gadfly, CSV
df = DataFrame(
    times=Float64[], accs=Int[], mean_accs=Float64[],
    median_accs=Float64[], max_equivs=String[], max_evals=Int[]
)


nruns = 10
for maxequiv in [1, 10, 50, 100, 500, 1000]
    for maxeval in [1, 5, 10, 25, 50]
        print("\rmaxPL=$(maxequiv), maxeval=$(maxeval)        ")

        opt_rt = @elapsed best_net, all_nets, _ = multi_search(deepcopy(gts[1]), q, 1; runs=nruns, maxeval=10000, maxequivPLs=maxequiv, opt_maxeval=maxeval)
        hwcds = [hardwiredClusterDistance(opt_net, net, false) for opt_net in all_nets]
        
        push!(df, [
            opt_rt, hardwiredClusterDistance(best_net, net, false),
            mean(hwcds), median(hwcds), "maxequivPLs=$(maxequiv)", maxeval
        ])
        CSV.write(joinpath(@__DIR__, "10run_param_perf.csv"), df)

        #p = Gadfly.plot(df, x=:max_equivs, y=:times, color=:max_evals, Geom.point)
        # stack_df = stack(df, [:times, :accs], variable_name="metric", value_name="value")
        # p = Gadfly.plot(stack_df, x=:max_evals, y=:value, xgroup=:max_equivs, ygroup=:metric, Geom.subplot_grid(Geom.boxplot))
        # display(p)
    end
end





