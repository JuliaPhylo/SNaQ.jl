using PhyloNetworks, Plots, SNaQ, DataFrames, PhyloCoalSimulations
include(joinpath(@__DIR__, "../../src/network_moves/rNNI_moves.jl"))
include(joinpath(@__DIR__, "../../src/network_moves/rSPR_moves.jl"))
include(joinpath(@__DIR__, "../../src/gradient_optimization/opt_API.jl"))
include(joinpath(@__DIR__, "../../src/gradient_optimization/search_API.jl"))


### Network example
net, gts = readnewick(joinpath(@__DIR__, "n1.netfile")), readmultinewick(joinpath(@__DIR__, "n1_10kgts.treefile")); net.isrooted = false;
q, t = countquartetsintrees(gts, showprogressbar=false);
df = readtableCF(DataFrame(tablequartetCF(q, t)));
snaq_net = SNaQ.readnewicklevel1(writenewick(net))
tre0 = readnewick(writenewick(gts[1]));


opt_rt = @elapsed opt_net, logPLs = search(tre0, q, net.numhybrids; seed=4, maxequivPLs = 100)
hardwiredClusterDistance(net, opt_net, false)





for j = 1:100
    rng = Random.seed!(j)

    net = generate_net(7, 1, j)
    gts = simulatecoalescent(net, 100, 1)
    semidirect_network!(net)
    q, t = countquartetsintrees(gts, showprogressbar=false)

    before_L = compute_loss(net, q)
    optimize_bls!(net, find_quartet_equations(net)[1], q, maxeval = 10)
    after_L = compute_loss(net, q)

    after_L > before_L || error("j = $(j)")
end

