"""
    fittedquartetCF(d::DataCF, format::Symbol)

Data frame with the observed and expected quartet concordance factors after
estimation of a network with `snaq!`, or fitting of quartet CF data on
a fixed network.
The format can be `:wide` (default) or `:long`.

- if `wide`, the output has one row per 4-taxon set, and each row has 10 columns: 4 columns
  for the taxon names, 3 columns for the observed CFs and 3 columns for the expected CF.
- if `long`, the output has one row per quartet, i.e. 3 rows per 4-taxon sets, and 7 columns:
  4 columns for the taxon names, one column to give the quartet resolution, one column for
  the observed CF and the last column for the expected CF.

See also: [`optimize!`](@ref) and [`loglik!`](@ref)
to update the fitted quartet CF expected
under a specific network, inside the DataCF object `d`.
"""
function fittedquartetCF(d::DataCF, format=:wide::Symbol)
    if format == :wide
        df=DataFrame(
                 tx1 = [q.taxon[1] for q in d.quartet],
                 tx2 = [q.taxon[2] for q in d.quartet],
                 tx3 = [q.taxon[3] for q in d.quartet],
                 tx4 = [q.taxon[4] for q in d.quartet],
                obsCF12=[q.obsCF[1] for q in d.quartet],
                obsCF13=[q.obsCF[2] for q in d.quartet],
                obsCF14=[q.obsCF[3] for q in d.quartet],
                expCF12=[q.expCF[1] for q in d.quartet],
                expCF13=[q.expCF[2] for q in d.quartet],
                expCF14=[q.expCF[3] for q in d.quartet])
    elseif format == :long
        nQ = length(d.quartet)
        df = DataFrame(
            tx1 = repeat([q.taxon[1] for q in d.quartet], inner=[3]),
            tx2 = repeat([q.taxon[2] for q in d.quartet], inner=[3]),
            tx3 = repeat([q.taxon[3] for q in d.quartet], inner=[3]),
            tx4 = repeat([q.taxon[4] for q in d.quartet], inner=[3]),
            quartet = repeat(["12_34","13_24","14_23"], outer=[nQ]),
            obsCF = Array{Float64}(undef, 3*nQ),
            expCF = Array{Float64}(undef, 3*nQ)  )
        row = 1
        for i in 1:nQ
            for j in 1:3
                df[row, 6] = d.quartet[i].obsCF[j]
                df[row, 7] = d.quartet[i].expCF[j]
                row += 1
            end
        end
    else
        error("format $(format) was not recognized. Should be :wide or :long.")
    end
    return df
end