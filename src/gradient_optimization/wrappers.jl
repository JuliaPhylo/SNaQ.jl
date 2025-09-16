

function snaq!(
    currT0::HybridNetwork,
    d::DataCF;
    hmax::Integer=1,
    liktolAbs::Float64=likAbs,
    Nfail::Integer=numFails,
    ftolRel::Float64=fRel,
    ftolAbs::Float64=fAbs,
    xtolRel::Float64=xRel,
    xtolAbs::Float64=xAbs,
    verbose::Bool=false,
    closeN::Bool=true,
    Nmov0::Vector{Int}=numMoves,
    runs::Integer=10,
    outgroup::AbstractString="none",
    filename::AbstractString="snaq",
    seed::Integer=0,
    probST::Float64=0.3,
    updateBL::Bool=true,
    probQR::Float64=0.0,
    qtolAbs::Float64=qAbs,
    qinfTest::Bool=false,
    propQuartets::Float64=1.0
)

    # kwargs not implemented that should be implemented:
    # - `filename`
    # - `updateBL`
    # - `probQR`
    # - `propQuartets`
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
        restrictions=restrictionset(max_level=1),
        ftolRel=ftolRel,
        ftolAbs=ftolAbs,
        xtolRel=xtolRel,
        xtolAbs=xtolAbs,
        propQuartets=propQuartets,
        filename=filename
    )[1]

end