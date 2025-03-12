using PhyloNetworks, Plots, SNaQ, DataFrames, PhyloCoalSimulations, Test, Random
include("../../src/network_moves/rNNI_moves.jl")
include("../../src/network_moves/rSPR_moves.jl")
include("../../src/gradient_optimization/opt_API.jl")
using PhyloPlots    # can be helpful, not always required


function get_data(L::Float64=0.5, ngt::Int=5_000, seed::Int=42)
    Random.seed!(seed)
    net = readnewick(joinpath(@__DIR__, "n1.netfile"))
    for E in net.edge E.length = (E.length > 0) ? L : 0 end
    net.isrooted = false
    rootonedge!(net, getroot(net).edge[1])
    gts = simulatecoalescent(net, ngt, 1);
    q, t = countquartetsintrees(gts, showprogressbar=false);
    semidirect_network!(net)
    return net, q
end



# Parameter 12 (the immediately parent to e and f) has a VERY large gradient - why??
net, q = get_data(1.0);

opt_net = deepcopy(net)
optimize_bls!(opt_net, q)
eqns, _, opt_params, idx_obj_map, _ = find_quartet_equations(opt_net);

@test sum((opt_params .- params).^2) < 0.1





