include("misc.jl")
include("CF_struct.jl")
const PN = PhyloNetworks;


function compute_logPL(net::HybridNetwork, obsCFs)
    quartet_eqns, _, _ = find_quartet_equations(net)
    return compute_logPL(quartet_eqns, net.edge, obsCFs)
end

function compute_logPL(quartet_eqns::Array{Tuple{Vector{eCFContribution},Vector{eCFContribution},Vector{eCFContribution}}}, edges::AbstractVector{PN.Edge}, oCFs)
    total_logPL = 0.0
    for (iter_eCF_eqns, iter_oCFs) in zip(quartet_eqns, oCFs)
        for (eqns, oCF) in zip(iter_eCF_eqns, iter_oCFs.data[1:3])
            if oCF == 0
                total_logPL += 0.0
            else
                eCF = compute_eCF(eqns, edges)
                total_logPL += 100 * oCF * log(eCF / oCF)
            end
        end
    end
    return total_logPL
end

function compute_eCFs(quartet_eqns::Array{Tuple{Vector{eCFContribution},Vector{eCFContribution},Vector{eCFContribution}}}, edges::AbstractVector{PN.Edge})
    eCFs = Array{Vector{Float64}}(undef, length(quartet_eqns))
    for (j, (eCF12_34_eqns, eCF13_24_eqns, eCF14_23_eqns)) in enumerate(quartet_eqns)
        eCFs[j] = [
            compute_eCF(eCF12_34_eqns, edges),
            compute_eCF(eCF13_24_eqns, edges),
            compute_eCF(eCF14_23_eqns, edges)
        ]
    end
    return eCFs
end

function compute_eCF(quartet_eqns::Vector{eCFContribution}, edges::AbstractVector{PN.Edge})
    total_eCF = 0.0
    for eqn in quartet_eqns
        total_eCF += compute_partial_eCF(eqn, edges)
    end
    return total_eCF
end

function compute_partial_eCF(eqn::eCFContribution, edges::AbstractVector{PN.Edge})
    exp_portion = (eqn.from_displayed_major) ? 1 - 2/3*exp(-sum(edges[e_idx].length for e_idx in eqn.internal_edges)) : 1/3*exp(-sum(edges[e_idx].length for e_idx in eqn.internal_edges))
    gamma_prod = (length(eqn.hyb_edges) == 0) ? 1 : prod(edges[e_idx].gamma for e_idx in eqn.hyb_edges)
    return gamma_prod * exp_portion
end

function compute_gradient(quartet_eqns::Array{Tuple{Vector{eCFContribution},Vector{eCFContribution},Vector{eCFContribution}}}, edges::AbstractVector{PN.Edge}, edge_and_hybrid_to_idx_map::Dict, oCFs)
    grad = fill(0.0, length(edge_and_hybrid_to_idx_map))
    for ((eqns1, eqns2, eqns3), iter_oCFs) in zip(quartet_eqns, oCFs)
        add_gradient_contributions!(eqns1, iter_oCFs.data[1], edges, edge_and_hybrid_to_idx_map, grad)
        add_gradient_contributions!(eqns2, iter_oCFs.data[2], edges, edge_and_hybrid_to_idx_map, grad)
        add_gradient_contributions!(eqns3, iter_oCFs.data[3], edges, edge_and_hybrid_to_idx_map, grad)
    end
    return grad
end

function add_gradient_contributions!(eqns::Vector{eCFContribution}, oCF, edges::AbstractVector{PN.Edge}, edge_and_hybrid_to_idx_map::Dict, grad::AbstractArray{Float64})
    
    eCF = compute_eCF(eqns, edges)
    
    # Compute gradient w.r.t. each branch length
    eCF_branch_prime = compute_eCF_branch_derivative(eqns, edges, oCF)
    edges_involved = union(reduce(vcat, eqn.internal_edges for eqn in eqns))
    for E_idx in edges_involved
        grad_idx = edge_and_hybrid_to_idx_map[edges[E_idx]]
        grad[grad_idx] += oCF / eCF * eCF_branch_prime
    end

    # Compute gradient w.r.t. gammas
    hybs_involved = union(reduce(vcat, eqn.hyb_edges for eqn in eqns))
    hybs_involved = [edge_idx for edge_idx in hybs_involved if edges[edge_idx].ismajor]
    for (j, H_idx) in enumerate(hybs_involved)
        eCF_gamma_prime = compute_eCF_gamma_derivative(eqns, edges, oCF, H_idx)
        grad_idx = edge_and_hybrid_to_idx_map[getchild(edges[H_idx])]
        grad[grad_idx] += oCF / eCF * eCF_gamma_prime
    end
end

function compute_eCF_gamma_derivative(eqns::Vector{eCFContribution}, edges::AbstractVector{PN.Edge}, oCF, H_idx::Int)
    total_deriv = 0.0
    for eqn in eqns
        partial_eCF = compute_partial_eCF(eqn, edges)
        eqn_has_edge = H_idx in eqn.hyb_edges
        denom = eqn_has_edge ? edges[H_idx].gamma : 1-edges[H_idx].gamma
        mult_neg_one = edges[H_idx].ismajor == eqn_has_edge
        total_deriv += partial_eCF / denom * (mult_neg_one ? -1 : 1)
    end
    return total_deriv
end

function compute_eCF_branch_derivative(eqns::Vector{eCFContribution}, edges::AbstractVector{PN.Edge}, oCF)
    total_deriv = 0.0
    for eqn in eqns
        expsumBL = exp(-sum(edges[E_idx].length for E_idx in eqn.internal_edges))
        gamma_prod = (length(eqn.hyb_edges) == 0) ? 1 : prod(edges[e_idx].gamma for e_idx in eqn.hyb_edges)
        exponential_portion = -1/3*expsumBL
        if eqn.from_displayed_major exponential_portion *= -2 end
        total_deriv += gamma_prod * exponential_portion
    end
    return total_deriv
end

function compute_eCF_branch_branch_2nd_derivative(eqns::Vector{eCFContribution}, edges::AbstractVector{PN.Edge}, oCF)
    
end

function compute_hessian!(H::AbstractMatrix, quartet_eqns::Array{Tuple{Vector{eCFContribution},Vector{eCFContribution},Vector{eCFContribution}}}, edges::AbstractVector{PN.Edge}, edge_and_hybrid_to_idx_map::Dict, oCFs, ngamma::Int)

    H[[(j, j) for j = 1:size(H)]] .= compute_gradient(quartet_eqns, edges, edge_and_hybrid_to_idx_map, oCFs)

end


function incr_bitvec!(bv::BitVector)
    carry = bv[1]
    bv[1] = !bv[1]
    carry_idx = 2

    while carry && carry_idx <= length(bv)
        carry = bv[carry_idx]
        bv[carry_idx] = !bv[carry_idx]
        carry_idx += 1
    end
end