using PhyloNetworks, SNaQ, PhyloCoalSimulations, Test, Random, StatsBase
import SNaQ: Node, semidirectnetwork!, findquartetequations, computeexpectedCFs


function get_data(L::Float64=0.5, seed::Int=42)
    Random.seed!(seed)
    net = readnewick(joinpath(@__DIR__, "n1.netfile"))
    semidirectnetwork!(net)
    for E in net.edge E.length = (E.length > 0) ? L : 0 end
    q = computeexpectedCFs(net)
    return net, q
end



@testset "optimize! gets close to truth" begin
    for L in [0.1, 0.25, 0.5, 1.0, 2.0, 5.0]
        for seed = 1:10
            net, q = get_data(L, seed);
            params = SNaQ.gatherparams(net);

            opt_net = deepcopy(net)
            for E in opt_net.edge
                if E.length != 0.0 E.length = rand()*L end
            end
            for H in opt_net.hybrid
                γ = rand()
                γ = max(γ, 1.0 - γ)
                getparentedge(H).gamma = γ
                getparentedgeminor(H).gamma = 1.0 - γ
            end

            firstoptL = optimize!(SNaQ.deepcopynetwork(opt_net), q; maxeval=100000)
            eqns = findquartetequations(opt_net)[1];
            secondoptL = optimize!(opt_net, eqns, q, Inf; maxeval=100000)
            opt_params = SNaQ.gatherparams(opt_net);
            @test firstoptL ≈ secondoptL atol=1e-5
            if !(sum(mean((opt_params .- params).^2)) < 1.0)
                @info L
                @info (sum(mean((opt_params .- params).^2)))
            end

            @test sum(mean(abs.(opt_params .- params))) < 1.0
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
        SNaQ.semidirectnetwork!(net)
        q, t = countquartetsintrees(gts, showprogressbar=false)
        qstat = Array{Float64}(undef, length(q), 3)
        for j in eachindex(q)
            for k = 1:3
                qstat[j, k] = q[j].data[k]
            end
        end
        q = qstat

        before_L = SNaQ.computeloss(net, q)
        SNaQ.optimize!(net, SNaQ.findquartetequations(net)[1], q, maxeval = 10)
        after_L = SNaQ.computeloss(net, q)

        @test after_L > before_L
        ntested += 1
    end
end

@testset "-logPL strictly improves with ρ != 0.0" begin
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
        SNaQ.semidirectnetwork!(net)
        q, t = countquartetsintrees(gts, showprogressbar=false)
        qstat = Array{Float64}(undef, length(q), 3)
        for j in eachindex(q)
            for k = 1:3
                qstat[j, k] = q[j].data[k]
            end
        end
        q = qstat

        ρ = rand() < 0.25 ? 1.0 : rand()
        before_L = SNaQ.computeloss(net, q, ρ)
        SNaQ.optimize!(net, SNaQ.findquartetequations(net)[1], q, ρ; maxeval = 10)
        after_L = SNaQ.computeloss(net, q, ρ)

        @test after_L > before_L
        ntested += 1
    end
end


