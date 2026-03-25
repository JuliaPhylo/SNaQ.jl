

function snaq!(
    currT0::HybridNetwork,
    d::DataCF;
    hmax::Integer=1,
    liktolAbs::Float64=1e-8,
    Nfail::Integer=3000,
    ftolRel::Float64=1e-8,
    ftolAbs::Float64=1e-8,
    xtolRel::Float64=1e-8,
    xtolAbs::Float64=1e-8,
    verbose::Bool=false,
    closeN::Bool=true,
    Nmov0::Vector{Int}=Int[],
    runs::Integer=10,
    outgroup::AbstractString="none",
    filename::AbstractString="snaq",
    seed::Integer=rand(Int),
    probST::Float64=0.3,
    updateBL::Bool=true,
    probQR::Float64=0.0,
    qtolAbs::Float64=1e-4,
    qinfTest::Bool=false,
    propQuartets::Float64=1.0,
    restrictions::Function=SNaQ.knownidentifiable,
    kwargs...
)

    # kwargs not implemented that should be implemented:
    # - `filename`
    # - `updateBL`
    # - `probQR`
    # - `qinfTest`
    # - `qtolAbs`

    return multisearch(
        currT0,
        d,
        hmax;
        runs=runs,
        maxequivPLs=Nfail,
        verbose=verbose,
        seed=seed,
        probST=probST,
        outgroup=outgroup,
        restrictions=restrictions,
        ftolRel=ftolRel,
        ftolAbs=ftolAbs,
        xtolRel=xtolRel,
        xtolAbs=xtolAbs,
        propQuartets=propQuartets,
        filename=filename
    )[1]

end