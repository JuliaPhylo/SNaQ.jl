

function snaq!(
    currT0::HybridNetwork,
    d::DataCF;
    hmax::Integer=1,
    liktolAbs::Float64=1e-10,
    Nfail::Integer=3000,
    ftolRel::Float64=1e-10,
    ftolAbs::Float64=1e-10,
    xtolRel::Float64=1e-10,
    xtolAbs::Float64=1e-10,
    verbose::Bool=false,
    closeN::Bool=true,
    Nmov0::Vector{Int}=Int[],
    runs::Integer=10,
    outgroup::AbstractString="none",
    filename::AbstractString="snaq",
    seed::Integer=0,
    probST::Float64=0.3,
    updateBL::Bool=true,
    probQR::Float64=0.0,
    qtolAbs::Float64=1e-4,
    qinfTest::Bool=false,
    propQuartets::Float64=1.0,
    restrictions::Function=SNaQ.knownidentifiable
)

    # kwargs not implemented that should be implemented:
    # - `filename`
    # - `updateBL`
    # - `probQR`
    # - `qinfTest`
    # - `qtolAbs`

    return multi_search(
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