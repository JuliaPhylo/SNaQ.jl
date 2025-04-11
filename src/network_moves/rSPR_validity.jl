

"""
Checks whether the rSPR move utilizing the given nodes is valid
according to the naming conventions in Figure 6 of [this paper](https://doi.org/10.1371/journal.pcbi.1005611).
"""
function is_valid_rSPR(w::Node, x::Node, y::Node, z::Node, xprime::Node, yprime::Node)
    # If type 1 (upper half of Fig 6): w must be child of z
    # If type 2 (lower ...): z must be child of w
    z.hybrid ? z in getchildren(w) : w in getchildren(z) || return false

    # y' must always be child of x'
    yprime in getchildren(xprime) || return false

    # all 6 nodes are unique
    if length(unique([w, x, y, z, xprime, yprime])) != 6
        (length(unique([w, x, y, z, xprime, yprime])) == 5 && (y == xprime || x == xprime)) || return false
    end

    if z.hybrid !(is_descendant_of(w, yprime)) || return false end
    if !z.hybrid && w.hybrid !(is_descendant_of(xprime, w)) || return false end

    # See Lemma 6, these conditions are necessary for the move to be valid
    if z.hybrid
        !is_descendant_of(w, yprime) || return false
    else
        !is_descendant_of(xprime, w) || return false
    end

    # Must be valid, then
    return true
end