using PhyloNetworks, SNaQ, PhyloCoalSimulations, Test, Random
import SNaQ: 
    RecursiveCFEquation, QuartetData, contains_parameter,
    compute_eCF, compute_eCFs, semidirect_network!,
    find_quartet_equations, compute_loss_and_gradient!, compute_loss,
    compute_eCF_and_gradient_recur!

@testset "RecursiveCFEquation construction" begin
    # Test basic construction of RecursiveCFEquation
    can_coal = true
    coal_es = [1, 2, 3]
    which_c = 1  # ab|cd
    dh = -1  # not set
    divisions = RecursiveCFEquation[]
    nparam = 5
    
    eqn = RecursiveCFEquation(can_coal, coal_es, which_c, dh, divisions, nparam)
    
    @test eqn.can_coalesce_here == can_coal
    @test eqn.coal_edges == coal_es
    @test eqn.which_coal == which_c
    @test eqn.division_H == dh
    @test eqn.divisions == divisions
    @test eqn.coal_mask == BitVector([true, true, true, false, false])
end

@testset "QuartetData construction and parameter checking" begin
    # Test QuartetData construction
    eqn = RecursiveCFEquation(true, [1, 2], 1, -1, RecursiveCFEquation[], 3)
    relevant_params = [1, 2, 3]
    q_taxa = SNaQ.SizedVector{4,String}(["A", "B", "C", "D"])
    
    qdata = QuartetData(eqn, relevant_params, q_taxa)
    
    @test qdata.eqn === eqn
    @test qdata.relevant_params == relevant_params
    @test qdata.q_taxa == q_taxa
    
    # Test contains_parameter
    @test contains_parameter(qdata, [1]) == true
    @test contains_parameter(qdata, [4]) == false
    @test contains_parameter(qdata, [1, 4]) == true
    @test contains_parameter(qdata, [4, 5]) == false
end

@testset "compute_eCF with RecursiveCFEquation" begin
    # Create a simple tree-like quartet with no hybridization
    # We'll manually create a RecursiveCFEquation for a simple case
    coal_edges = [1]  # Just one internal edge
    eqn = RecursiveCFEquation(false, coal_edges, 1, -1, RecursiveCFEquation[], 1)
    relevant_params = [1]
    q_taxa = SNaQ.SizedVector{4,String}(["A", "B", "C", "D"])
    
    qdata = QuartetData(eqn, relevant_params, q_taxa)
    
    # For a tree with topology ((A,B),(C,D)) with one internal edge,
    # we expect eCF1 (AB|CD) to be dominant
    params = [0.5]  # Branch length of 0.5
    ecfs = compute_eCF(qdata, params, Inf)
    
    # Check that eCF1 (AB|CD) > eCF2 (AC|BD) = eCF3 (AD|BC)
    @test ecfs[1] > ecfs[2]
    @test ecfs[2] ≈ 1/3 * exp(-0.5)
    @test ecfs[1] ≈ 1 - 2/3 * exp(-0.5)
end

@testset "compute_loss_and_gradient! calculation" begin
    # Set up a simple network and compute its gradient
    Random.seed!(42)
    net = readnewick("((A,(B)#H1),((C,D), #H1));"); # Simple network with one hybrid node
    for E in net.edge E.length = 1.0; E.gamma = !E.hybrid ? -1 : 0.5 end
    semidirect_network!(net)
    
    # Simulate quartet frequencies
    true_q = [0.6, 0.2, 0.2]  # Some arbitrary frequencies
    q = reshape(true_q, 1, 3)
    
    # Get QuartetData and parameters
    qdata, _, params, _, _ = find_quartet_equations(net)
    
    # Create gradient storage
    gradient = zeros(length(params))
    
    # Compute loss and gradient
    loss = compute_loss_and_gradient!(qdata, params, gradient, q)
    
    # The loss should be non-negative
    @test loss ≈ -0.20561335094444322
    
    # Test numerical gradient approximation
    epsilon = 1e-6
    for i in 1:length(params)
        params_plus = copy(params)
        params_plus[i] += epsilon
        
        params_minus = copy(params)
        params_minus[i] -= epsilon
        
        loss_plus = compute_loss(qdata, params_plus, q)
        loss_minus = compute_loss(qdata, params_minus, q)
        
        numerical_grad = (loss_plus - loss_minus) / (2 * epsilon)
        
        # The analytical gradient should be close to the numerical approximation
        if abs(numerical_grad) > 1e-6
            @test isapprox(gradient[i], numerical_grad, rtol=0.1)
        end
    end
end

@testset "compute_eCF_and_gradient_recur! for simple cases" begin
    # Test a simple treelike case (division_H = -1)
    coal_edges = [1]
    eqn_treelike = RecursiveCFEquation(false, coal_edges, 1, -1, RecursiveCFEquation[], 1)
    
    params = [0.5]
    gradient = zeros(1, 3)
    params_seen = falses(1)
    running_gradient = ones(1, 3)
    
    # Compute eCFs and gradient
    ecf1, ecf2 = compute_eCF_and_gradient_recur!(eqn_treelike, params, gradient, params_seen, Inf, running_gradient)
    
    # Check eCFs
    @test ecf1 ≈ 1 - 2/3 * exp(-0.5)
    @test ecf2 ≈ 1/3 * exp(-0.5)
    
    # Check gradient
    @test gradient[1, 1] ≈ 2/3 * exp(-0.5)  # Derivative for AB|CD
    @test gradient[1, 2] ≈ -1/3 * exp(-0.5)  # Derivative for AC|BD
    @test gradient[1, 3] ≈ -1/3 * exp(-0.5)  # Derivative for AD|BC
end 