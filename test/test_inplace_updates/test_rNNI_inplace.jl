using PhyloNetworks
using Random, PhyloCoalSimulations
using Test
using SNaQ
import SNaQ:
    QuartetData, semidirectnetwork!, findquartetequations, gatheroptimizationinfo,
    allvalidrNNI1nodes, performrNNI1!, applyrNNI1update!,
    allvalidrNNI2nodes, performrNNI2!, applyrNNI2update!, sameparameterlayout
include("misc.jl")




########### GENERATE DATA ###########
Random.seed!(42)
net = generate_net(15, 2);
rootonedge!(net, getroot(net).edge[1]);
gts = simulatecoalescent(net, 100, 1);
q = SNaQ.gatherCFmatrix(readtrees2CF(gts; writeTab=false, writeSummary=false));
semidirectnetwork!(net);
net_qdata, net_param_map, net_params, _ = findquartetequations(net);

outofplace_times = [];
inplace_times = [];
outofplace_losses = [];
inplace_losses = [];
n_moves = length(allvalidrNNI1nodes(net));
for j = 1:n_moves
    print("\r$(j)/$(n_moves)")
    # moves that change which edges carry parameters are never updated in place by `search`
    net_copy = deepcopy(net); performrNNI1!(net_copy, allvalidrNNI1nodes(net_copy)[j]...)
    sameparameterlayout(net, net_copy) || continue
    # in-place
    inplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = allvalidrNNI1nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        performrNNI1!(net_copy, s, t, u, v);
        _, copy_param_map, _, copy_params, _ = gatheroptimizationinfo(net_copy);

        copy_qdata = Array{QuartetData}(undef, length(net_qdata));
        applyrNNI1update!(net_copy, net_qdata, copy_qdata, copy_param_map, u, v, 0.0)
        push!(inplace_losses, computeSNaQscore!(copy_qdata, copy_params, q))
    end
    push!(inplace_times, inplace_time)


    # not in-place
    outofplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = allvalidrNNI1nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        performrNNI1!(net_copy, s, t, u, v);
        correct_qdata, correct_param_map, correct_params, _ = findquartetequations(net_copy)
        push!(outofplace_losses, computeSNaQscore!(correct_qdata, correct_params, q))
    end
    push!(outofplace_times, outofplace_time)
end
@test all(isapprox.(outofplace_losses, inplace_losses))




outofplace_times = [];
inplace_times = [];
outofplace_losses = [];
inplace_losses = [];
n_moves = length(allvalidrNNI2nodes(net));
for j = 1:n_moves
    print("\r$(j)/$(n_moves)")
    # moves that change which edges carry parameters are never updated in place by `search`
    net_copy = deepcopy(net); performrNNI2!(net_copy, allvalidrNNI2nodes(net_copy)[j]...)
    sameparameterlayout(net, net_copy) || continue
    # in-place
    inplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = allvalidrNNI2nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        performrNNI2!(net_copy, s, t, u, v);
        _, copy_param_map, _, copy_params, _ = gatheroptimizationinfo(net_copy);

        copy_qdata = Array{QuartetData}(undef, length(net_qdata));
        applyrNNI2update!(net_copy, net_qdata, copy_qdata, copy_param_map, s, t, u, v, 0.0)
        push!(inplace_losses, computeSNaQscore!(copy_qdata, copy_params, q))
    end
    push!(inplace_times, inplace_time)


    # not in-place
    outofplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = allvalidrNNI2nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        performrNNI2!(net_copy, s, t, u, v);
        correct_qdata, correct_param_map, correct_params, _ = findquartetequations(net_copy)
        push!(outofplace_losses, computeSNaQscore!(correct_qdata, correct_params, q))
    end
    push!(outofplace_times, outofplace_time)
end
@test all(isapprox.(outofplace_losses, inplace_losses))


