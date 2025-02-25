using PhyloNetworks



function is_valid_rSPR(w::Node, x::Node, y::Node, z::Node, xprime::Node, yprime::Node)
    # Type 1
    z in getchildren(x) &&
        y in getchildren(z) &&
        w in getchildren(z) &&
        yprime in getchildren(xprime) &&
        !x.hybrid && !z.hybrid &&
        (!xprime.hybrid || xprime.hybrid != yprime.hybrid) &&
        return true

    # Type 2
    z in getchildren(x) &&
        z in getchildren(w) &&
        y in getchildren(z) &&
        yprime in getchildren(xprime) &&
        z.hybrid && !x.hybrid && !w.hybrid && !y.hybrid &&
        (!xprime.hybrid || xprime.hybrid != yprime.hybrid) &&
        return true

    # Invalid
    return false
end