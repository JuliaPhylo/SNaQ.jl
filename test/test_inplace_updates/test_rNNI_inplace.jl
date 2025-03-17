using PhyloNetworks
using Random, PhyloCoalSimulations
using Test
include("../../src/gradient_optimization/opt_API.jl")
include("../../src/gradient_optimization/search_API.jl")
include("../../src/gradient_optimization/CF_recursive_blocks.jl")
include("../../src/gradient_optimization/inplace_updates.jl")
include("misc.jl")




########### GENERATE DATA ###########
Random.seed!(42)
net = generate_net(15, 2);
rootonedge!(net, getroot(net).edge[1]);
gts = simulatecoalescent(net, 100, 1);
q, t = countquartetsintrees(gts, showprogressbar=false);
semidirect_network!(net);
net_qdata, net_param_map, net_params, _ = find_quartet_equations(net);

outofplace_times = [];
inplace_times = [];
outofplace_losses = [];
inplace_losses = [];
n_moves = length(all_valid_rNNI1_nodes(net));
for j = 1:n_moves
    print("\r$(j)/$(n_moves)")
    # in-place
    inplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = all_valid_rNNI1_nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        perform_rNNI1!(net_copy, s, t, u, v);
        _, _, _, copy_params, _ = gather_optimization_info(net_copy);

        copy_qdata = Array{QuartetData}(undef, length(net_qdata));
        apply_rNNI1_update!(net_copy, net_qdata, copy_qdata, net_param_map, u, Inf)
        push!(inplace_losses, compute_loss(copy_qdata, copy_params, q))
    end
    push!(inplace_times, inplace_time)


    # not in-place
    outofplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = all_valid_rNNI1_nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        perform_rNNI1!(net_copy, s, t, u, v);
        correct_qdata, correct_param_map, correct_params, _ = find_quartet_equations(net_copy)
        push!(outofplace_losses, compute_loss(correct_qdata, correct_params, q))
    end
    push!(outofplace_times, outofplace_time)
end
@test all(outofplace_losses .== inplace_losses)




outofplace_times = [];
inplace_times = [];
outofplace_losses = [];
inplace_losses = [];
n_moves = length(all_valid_rNNI2_nodes(net));
for j = 1:n_moves
    print("\r$(j)/$(n_moves)")
    # in-place
    inplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = all_valid_rNNI2_nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        perform_rNNI2!(net_copy, s, t, u, v);
        _, _, _, copy_params, _ = gather_optimization_info(net_copy);

        copy_qdata = Array{QuartetData}(undef, length(net_qdata));
        apply_rNNI2_update!(net_copy, net_qdata, copy_qdata, net_param_map, s, t, u, v, Inf)
        push!(inplace_losses, compute_loss(copy_qdata, copy_params, q))
    end
    push!(inplace_times, inplace_time)


    # not in-place
    outofplace_time = @elapsed begin
        net_copy = deepcopy(net);
        Random.seed!(42)
        valid_stuvs = all_valid_rNNI2_nodes(net_copy);
        s, t, u, v = valid_stuvs[j];
        perform_rNNI2!(net_copy, s, t, u, v);
        correct_qdata, correct_param_map, correct_params, _ = find_quartet_equations(net_copy)
        push!(outofplace_losses, compute_loss(correct_qdata, correct_params, q))
    end
    push!(outofplace_times, outofplace_time)
end
@test all(outofplace_losses .== inplace_losses)


