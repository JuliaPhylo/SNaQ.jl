using PhyloNetworks, SNaQ, PhyloCoalSimulations, Test, Random, StatsBase
import SNaQ: Node, semidirect_network!, find_quartet_equations, compute_eCFs


function get_data(L::Float64=0.5, seed::Int=42)
    Random.seed!(seed)
    net = readnewick(joinpath(@__DIR__, "n1.netfile"))
    semidirect_network!(net)
    for E in net.edge E.length = (E.length > 0) ? L : 0 end
    q = compute_eCFs(net)
    return net, q
end



@testset "optimize_bls! gets close to truth" begin
    for L in [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
        for seed = 1:10
            net, q = get_data(L, seed);
            _, _, params, _, _ = find_quartet_equations(net);

            opt_net = deepcopy(net)
            for E in opt_net.edge
                if E.length != 0.0 E.length = rand()*L end
                if E.gamma != -1.0 E.gamma = 0.5 end
            end

            optimize_bls!(opt_net, q; maxeval=100000)
            eqns, _, opt_params, idx_obj_map, _ = find_quartet_equations(opt_net);
            if !(sum(mean((opt_params .- params).^2)) < 1.0)
                @info L
                @info (sum(mean((opt_params .- params).^2)))
            end

            @test sum(mean((opt_params .- params).^2)) < 1.0
        end
    end
end

@testset "-logPL strictly improves" begin
    ntested::Int = 0
    nfail::Int = 0
    while ntested < 100
        rng = Random.seed!(ntested + nfail)

        net = generate_net(7, 1, ntested + nfail)
        while shrink2cycles!(net) || shrink3cycles!(net) continue end
        if net.numhybrids == 0
            nfail += 1
            continue
        end
        gts = simulatecoalescent(net, 100, 1)
        SNaQ.semidirect_network!(net)
        q, t = countquartetsintrees(gts, showprogressbar=false)
        qstat = Array{Float64}(undef, length(q), 3)
        for j in eachindex(q)
            for k = 1:3
                qstat[j, k] = q[j].data[k]
            end
        end
        q = qstat

        before_L = SNaQ.compute_loss(net, q)
        SNaQ.optimize_bls!(net, SNaQ.find_quartet_equations(net)[1], q, maxeval = 10)
        after_L = SNaQ.compute_loss(net, q)

        @test after_L > before_L
        ntested += 1
    end
end


