using PhyloNetworks, SNaQ, DataFrames, PhyloCoalSimulations
using Test, Random

include(joinpath(@__DIR__, "../test_inplace_updates/misc.jl"))

### Network example
net = readnewick(joinpath(@__DIR__, "n1.netfile"))
q = SNaQ.compute_eCFs(net);
tre0 = majortree(net);


opt_rt = @elapsed opt_net, logPLs = search(tre0, q, net.numhybrids; seed=5, maxequivPLs = 1000)
@test hardwiredClusterDistance(net, opt_net, false) == 0




# error("@btime")
# using BenchmarkTools
# @btime search(tre0, q, net.numhybrids; seed=5, maxequivPLs = 1000)
# # 4.7s, 5.29 GiB
