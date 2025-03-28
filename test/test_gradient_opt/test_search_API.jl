using Revise


using PhyloNetworks, Plots, SNaQ, DataFrames
include("../../src/network_moves/rNNI_moves.jl")
include("../../src/network_moves/rSPR_moves.jl")
include("../../src/gradient_optimization/opt_API.jl")
include("../../src/gradient_optimization/search_API.jl")
ENV["JULIA_DEBUG"] = Main
ENV["JULIA_DEBUG"] = nothing

using PhyloCoalSimulations, PhyloPlots


### Network example
net, gts = readnewick(joinpath(@__DIR__, "n1.netfile")), readmultinewick(joinpath(@__DIR__, "n1_10kgts.treefile")); net.isrooted = false;
q, t = countquartetsintrees(gts, showprogressbar=false);
df = readtableCF(DataFrame(tablequartetCF(q, t)));
snaq_net = SNaQ.readnewicklevel1(writenewick(net))
tre0 = readnewick(writenewick(gts[1]));
perform_rNNI1!(tre0, sample_rNNI_parameters(tre0, 1, Random.seed!(0))...);
perform_rNNI1!(tre0, sample_rNNI_parameters(tre0, 1, Random.seed!(0))...);
perform_rNNI1!(tre0, sample_rNNI_parameters(tre0, 1, Random.seed!(0))...);
perform_rNNI1!(tre0, sample_rNNI_parameters(tre0, 1, Random.seed!(0))...);


opt_rt = @elapsed opt_net, logPLs = search(tre0, q, net.numhybrids; seed=4, maxequivPLs = 100)
hardwiredClusterDistance(net, opt_net, false)







nprocs() < 5 && addprocs(5 - nprocs())
@everywhere include("../../src/gradient_optimization/search_API.jl")


function split_search(nsplit, maxiter, seed, t0, q, hmax; kwargs...)

    rng = Random.seed!(seed)
    best_net = deepcopy_network(t0)
    prev_PL = compute_loss(best_net, q)
    for iter = 1:maxiter
        seeds = rand(rng, Int, nsplit)
        res = Distributed.pmap(1:nsplit) do j
            opt_net, opt_PLs = search(best_net, q, net.numhybrids; seed=seeds[j], probST = 1.0)
            return opt_net, opt_PLs
        end
        best_net = res[findmax(np -> np[2], res)[2]][1]
        best_PL = findmax(np -> np[2], res)[1]

        if prev_PL > best_PL
            break
        end
        prev_PL = best_PL
    end
    return best_net, compute_loss(best_net, q)

end




split_hwcds = [];
split_rts = [];
multi_hwcds = [];
multi_rts = [];
while true
    print("j=$(length(split_hwcds)+1):\t")
    j = sample(1:length(gts))
    srt = @elapsed snet, sPL = split_search(nprocs(), 25, 100+length(split_hwcds), deepcopy(gts[j]), q, net.numhybrids; maxequivPLs=75, opt_maxeval=10)
    push!(split_rts, srt)
    push!(split_hwcds, hardwiredclusterdistance(snet, net, false))
    print("split=($(round(srt, digits=2)), $(split_hwcds[length(split_hwcds)]))")

    mrt = @elapsed mnet, _ = multi_search(deepcopy(gts[j]), q, net.numhybrids; seed=100+length(split_hwcds), runs=nprocs(), maxequivPLs=1000, opt_maxeval=10)
    push!(multi_rts, mrt)
    push!(multi_hwcds, hardwiredclusterdistance(mnet, net, false))

    println(", multi=($(round(mrt, digits=2)), $(multi_hwcds[length(multi_hwcds)]))")
end

snet, sPL = split_search(5, 25, 1, tre0, q, net.numhybrids; maxequivPLs = 50)
hardwiredclusterdistance(snet, net, false)

