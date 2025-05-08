using PhyloNetworks, SNaQ, PhyloCoalSimulations, Test, Random, StatsBase
import SNaQ: Node, semidirect_network!, find_quartet_equations, compute_eCFs


function get_data(L::Float64=0.5, seed::Int=42)
    Random.seed!(seed)
    net = readnewick(joinpath(@__DIR__, "n1.netfile"))
    for E in net.edge E.length = (E.length > 0) ? L : 0 end
    semidirect_network!(net)
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
            if !(sum(mean((opt_params .- params).^2)) < 0.5)
                @info L
            end

            @test sum(mean((opt_params .- params).^2)) < 0.5
        end
    end
end


