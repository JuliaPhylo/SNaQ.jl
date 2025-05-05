using PhyloNetworks, SNaQ, DataFrames, PhyloCoalSimulations
using Test, Random

include(joinpath(@__DIR__, "../test_inplace_updates/mics.jl"))

### Network example
net, gts = readnewick(joinpath(@__DIR__, "n1.netfile")), readmultinewick(joinpath(@__DIR__, "n1_10kgts.treefile")); net.isrooted = false;
q, t = countquartetsintrees(gts, showprogressbar=false);
df = readtableCF(DataFrame(tablequartetCF(q, t)));
snaq_net = SNaQ.readnewicklevel1(writenewick(net))
tre0 = readnewick(writenewick(gts[1]));


opt_rt = @elapsed opt_net, logPLs = search(tre0, q, net.numhybrids; seed=4, maxequivPLs = 100)
@test hardwiredClusterDistance(net, opt_net, false) == 0


@testset "multi-opt" begin
    for j = 1:100
        rng = Random.seed!(j)

        net = generate_net(7, 1, j)
        gts = simulatecoalescent(net, 100, 1)
        SNaQ.semidirect_network!(net)
        q, t = countquartetsintrees(gts, showprogressbar=false)
        qstat = Array{Float64}(undef, length(q), 3)
        for j = 1:length(q)
            for k = 1:3
                qstat[j, k] = q[j].data[k]
            end
        end
        q = qstat

        before_L = SNaQ.compute_loss(net, qstat)
        SNaQ.optimize_bls!(net, SNaQ.find_quartet_equations(net)[1], q, maxeval = 10)
        after_L = SNaQ.compute_loss(net, qstat)

        after_L > before_L || error("j = $(j)")
        @test after_L > before_L
    end
end

