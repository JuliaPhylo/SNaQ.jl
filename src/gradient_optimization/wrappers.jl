

# function snaq!(
#     currT0::HybridNetwork,
#     d::DataCF;
#     hmax::Integer=1,
#     liktolAbs::Float64=likAbs,
#     Nfail::Integer=numFails,
#     ftolRel::Float64=fRel,
#     ftolAbs::Float64=fAbs,
#     xtolRel::Float64=xRel,
#     xtolAbs::Float64=xAbs,
#     verbose::Bool=false,
#     closeN::Bool=true,
#     Nmov0::Vector{Int}=Int[],
#     runs::Integer=10,
#     outgroup::AbstractString="none",
#     filename::AbstractString="snaq",
#     seed::Integer=0,
#     probST::Float64=0.3,
#     updateBL::Bool=true,

#     # New arguments
#     restrictions::Function=restriction_set(;rooted_tree_child=true, galled_network=true)
# )

#     multi_search(
#         currT0,
#         d,        hmax;
#         runs=runs,
#         seed=seed,
#         logprefix=filename,
#         restrictions=restrictions,
#         probST=probST,
#     )

# end