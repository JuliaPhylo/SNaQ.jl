using PhyloNetworks, PhyloCoalSimulations
include("../../src/network_moves/rNNI_moves.jl")
include("../../src/network_moves/rSPR_moves.jl")
include("../../src/gradient_optimization/opt_API.jl")


function load_net_and_gts(; ngt=1000)
    net = readnewick("((a,b)i1,(((c1,c2),#H1)i6,((d1,d2),((e,f)i5)#H1)i3)i2);")
    for E in net.edge
        E.length = 0.5
        if E.hybrid
            E.gamma = 0.5
            if !E.ismajor E.length = 0.0 end
        end
    end
    gts = simulatecoalescent(net, ngt, 1);
    semidirect_network!(net)
    net.hybrid[1].name = "i4"   # do it here so that PhyloNetworks doesn't through a warning
    return net, gts
end

net, gts = load_net_and_gts(; ngt=10000);
q, t = countquartetsintrees(gts, showprogressbar=false);
net0 = deepcopy(net)
optimize_bls!(net0, q)

PhyloPlots.plot(net, shownodelabel=true, showedgelength=true)
PhyloPlots.plot(net0, shownodelabel=true, showedgelength=true)



