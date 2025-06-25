using SNaQ, PhyloNetworks, Random, Test, StatsBase, PhyloCoalSimulations
include("test_inplace_updates/misc.jl")
import SNaQ: Node, semidirect_network!, find_quartet_equations, compute_eCFs,
    is_valid_add_hybrid,
    add_hybrid!, remove_hybrid!,
    getparentedge, getparent,
    sample_move_reticulate_origin_parameters,
    is_valid_move_reticulate_origin, move_reticulate_origin!,deepcopy_network, semidirect_network!,
    generate_move_proposal, apply_move!,
    move_reticulate_target!, move_reticulate_origin!,
    semidirect_network!,
    perform_rNNI!, perform_rNNI1!, perform_rNNI2!, perform_rNNI3!, perform_rNNI4!,
    is_valid_rNNI1, is_valid_rNNI2, is_valid_rNNI3, is_valid_rNNI4,
    all_valid_rNNI_nodes, perform_random_rNNI!, semidirect_network!,
    all_valid_rNNI1_nodes, all_valid_rNNI2_nodes,
    apply_move!, perform_rSPR!, is_valid_rSPR, semidirect_network!,
    sample_rSPR_parameters

test_files = [
    "test_gradient_opt/test_opt_API.jl",
    "test_gradient_opt/test_search_API.jl",
    "test_gradient_opt/test_CF_recursive_blocks.jl",
    "test_network_moves/test_add_remove_retic.jl",
    "test_network_moves/test_move_target_origin.jl",
    "test_network_moves/test_rNNI_moves.jl",
    "test_network_moves/test_rSPR_moves.jl",
    "test_network_moves/test_all_moves.jl"
]

function redirect_outputs(f)
    redirect_stdout(devnull) do
        f()
    end
end

@testset "SNaQ.jl" begin
    for file in test_files
        @testset "$file" begin
            redirect_outputs() do
                include(file)
            end
        end
    end
end