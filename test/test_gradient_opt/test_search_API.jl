using PhyloNetworks, Plots, SNaQ, DataFrames
include("../../src/network_moves/rNNI_moves.jl")
include("../../src/network_moves/rSPR_moves.jl")
include("../../src/gradient_optimization/opt_API.jl")
include("../../src/gradient_optimization/search_API.jl")

### Network example
net, gts = readnewick(joinpath(@__DIR__, "n1.netfile")), readmultinewick(joinpath(@__DIR__, "n1_10kgts.treefile")); net.isrooted = false;
q, t = countquartetsintrees(gts, showprogressbar=false);
df = readtableCF(DataFrame(tablequartetCF(q, t)));
tre0 = readnewick(writenewick(gts[1]));
perform_random_rNNI!(tre0);

opt_rt = @elapsed opt_net, logPLs = search(tre0, q, net.numhybrids; maxeval = 1000, maxequivPLs = 100)
hardwiredClusterDistance(net, opt_net, false)

snaq_rt = @elapsed snaq_net = snaq!(tre0, df, hmax=net.numhybrids, runs=1);
hardwiredClusterDistance(net, snaq_net, false)


compute_logPL(opt_net, q), compute_logPL(net, q)

Plots.plot(1:length(logPLs), logPLs)

PhyloPlots.plot(net);
PhyloPlots.plot(opt_net);



### Tree example
tre = readnewick(joinpath(@__DIR__, "n1.netfile"));
gts = simulatecoalescent(tre, 1000, 1);
tre.isrooted = false;
q, _ = countquartetsintrees(gts, showprogressbar=false);
tre0 = readnewick(writenewick(tre));
tre0.isrooted = false;
perform_random_rNNI!(tre0);perform_random_rNNI!(tre0);perform_random_rNNI!(tre0);

opt_tre, logPLs = search(tre0, q, tre.numhybrids; maxeval = 1000, maxequivPLs = 100)
hardwiredClusterDistance(tre, opt_tre, false)

PhyloPlots.plot(tre)
PhyloPlots.plot(tre0)


### hmax = 3 will probably give us plenty of fun, new errors to debug...
tre = readnewick(joinpath(@__DIR__, "n1.netfile"));
gts = simulatecoalescent(tre, 1000, 1);
tre.isrooted = false;
q, _ = countquartetsintrees(gts, showprogressbar=false);
tre0 = readnewick(writenewick(tre));
tre0.isrooted = false;
perform_random_rNNI!(tre0);perform_random_rNNI!(tre0);perform_random_rNNI!(tre0);

opt_tre, logPLs = search(tre0, q, 3; maxeval = 1000, maxequivPLs = 100)
hardwiredClusterDistance(tre, opt_tre, false)