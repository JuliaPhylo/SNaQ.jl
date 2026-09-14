# The running gradient carried down an equation tree.

"""
    RunningGradient(nparam)

Scratch space for [`computeexpectedCFandgradientrecur!`](@ref), holding the *running
gradient*: the product of factors accumulated from the root of an equation tree down to the
node being evaluated, for each parameter and each of the 3 quartet topologies.

Each step of the recursion scales the running gradient of every parameter by the same
factor, except for the γ of the reticulation being divided on and the branches the current
equation coalesces over. It is therefore stored as one shared triple (`shared`) plus the
handful of *exceptions* that have drifted from it, which keeps every sweep proportional to
the number of exceptions rather than to `nparam`.

One is built per thread by [`threadscratch`](@ref) and reused across quartets.
"""
mutable struct RunningGradient
    shared::Vector{Float64}         # length 3: running gradient of every non-exceptional parameter
    exceptionparam::Vector{Int}     # parameter index held in each exception slot
    exceptionvalue::Matrix{Float64} # 3 x capacity: running gradient of each exception
    exceptionslot::Vector{Int}    # nparam-long: parameter -> its slot, 0 if not an exception
    params_seen::Vector{Bool}       # nparam-long: has an ancestor already fixed this parameter?
    nexceptions::Int
    savestack::Vector{Float64}      # checkpoint arena for saverunninggradient!
    savetop::Int

    function RunningGradient(nparam::Int)
        # Both the exception list and the checkpoint arena grow on demand.
        cap = 8
        new(ones(3), zeros(Int, cap), Matrix{Float64}(undef, 3, cap),
            zeros(Int, nparam), zeros(Bool, nparam), 0, Vector{Float64}(undef, 64), 0)
    end
end


"""
Resets `rg` to the state a recursion starts from: running gradient 1 everywhere, no
exceptions, nothing seen.
"""
@inline function resetrunninggradient!(rg::RunningGradient)
    rg.shared[1] = 1.0; rg.shared[2] = 1.0; rg.shared[3] = 1.0
    # Every flag the recursion sets is on a parameter it also made an exception, so clearing
    # the exceptions clears `params_seen` too.
    @inbounds for s = 1:rg.nexceptions
        rg.exceptionslot[rg.exceptionparam[s]] = 0
        rg.params_seen[rg.exceptionparam[s]] = false
    end
    rg.nexceptions = 0
    rg.savetop = 0
    return nothing
end


"""
Running gradient of parameter `p` for topology `k`: its own value if `p` is an exception,
otherwise the shared value.
"""
@inline function runninggradient(rg::RunningGradient, p::Int, k::Int)::Float64
    @inbounds s = rg.exceptionslot[p]
    return s == Int(0) ? (@inbounds rg.shared[k]) : (@inbounds rg.exceptionvalue[k, s])
end


"""
Makes `p` an exception, starting from the value it currently has. A no-op if `p` already is one.
"""
@inline function addexception!(rg::RunningGradient, p::Int)
    @inbounds rg.exceptionslot[p] == Int(0) || return nothing
    if rg.nexceptions == length(rg.exceptionparam)
        newcap = 2 * rg.nexceptions
        resize!(rg.exceptionparam, newcap)
        neweval = Matrix{Float64}(undef, 3, newcap)
        @inbounds copyto!(view(neweval, :, 1:rg.nexceptions), rg.exceptionvalue)
        rg.exceptionvalue = neweval
    end
    rg.nexceptions += 1
    s = rg.nexceptions
    @inbounds begin
        rg.exceptionparam[s] = p
        rg.exceptionslot[p] = Int(s)
        rg.exceptionvalue[1, s] = rg.shared[1]; rg.exceptionvalue[2, s] = rg.shared[2]; rg.exceptionvalue[3, s] = rg.shared[3]
    end
    return nothing
end


"""
Drops every exception added since `rg.nexceptions` was `nex0`. Exceptions are added and
dropped in LIFO order, so this restores exactly the set that was live at that point.
"""
@inline function dropexceptionsto!(rg::RunningGradient, nex0::Int)
    @inbounds for s = rg.nexceptions:-1:(nex0 + 1)
        rg.exceptionslot[rg.exceptionparam[s]] = 0
    end
    rg.nexceptions = nex0
    return nothing
end


"""
Checkpoints the running gradient onto `rg.savestack`, returning a handle for
[`restorerunninggradient!`](@ref). Needed because a division's factor can be exactly 0, so
the scaling cannot simply be divided back out.
"""
@inline function saverunninggradient!(rg::RunningGradient)
    need = rg.savetop + 3 + 3 * rg.nexceptions
    length(rg.savestack) < need && resize!(rg.savestack, max(need, 2 * length(rg.savestack)))
    @inbounds begin
        rg.savestack[rg.savetop + 1] = rg.shared[1]
        rg.savestack[rg.savetop + 2] = rg.shared[2]
        rg.savestack[rg.savetop + 3] = rg.shared[3]
        o = rg.savetop + 3
        for s = 1:rg.nexceptions
            rg.savestack[o + 1] = rg.exceptionvalue[1, s]
            rg.savestack[o + 2] = rg.exceptionvalue[2, s]
            rg.savestack[o + 3] = rg.exceptionvalue[3, s]
            o += 3
        end
    end
    rg.savetop = need
    return need
end


"""
Restores the checkpoint [`saverunninggradient!`](@ref) returned as `top`.
"""
@inline function restorerunninggradient!(rg::RunningGradient, top::Int)
    base = top - 3 - 3 * rg.nexceptions
    @inbounds begin
        rg.shared[1] = rg.savestack[base + 1]
        rg.shared[2] = rg.savestack[base + 2]
        rg.shared[3] = rg.savestack[base + 3]
        o = base + 3
        for s = 1:rg.nexceptions
            rg.exceptionvalue[1, s] = rg.savestack[o + 1]
            rg.exceptionvalue[2, s] = rg.savestack[o + 2]
            rg.exceptionvalue[3, s] = rg.savestack[o + 3]
            o += 3
        end
    end
    rg.savetop = base
    return nothing
end


"""
Scales the running gradient of every parameter by `f`, except parameter `skip` (0 scales
everything).
"""
@inline function scalerunninggradient!(rg::RunningGradient, f::Float64, skip::Int=0)
    @inbounds begin
        rg.shared[1] *= f; rg.shared[2] *= f; rg.shared[3] *= f
        for s = 1:rg.nexceptions
            rg.exceptionparam[s] == skip && continue
            rg.exceptionvalue[1, s] *= f; rg.exceptionvalue[2, s] *= f; rg.exceptionvalue[3, s] *= f
        end
    end
    return nothing
end


"""
Scales the running gradient of the single exceptional parameter `p` by `f`.
"""
@inline function scaleexception!(rg::RunningGradient, p::Int, f::Float64)
    @inbounds s = rg.exceptionslot[p]
    s == Int(0) && return nothing
    @inbounds begin
        rg.exceptionvalue[1, s] *= f; rg.exceptionvalue[2, s] *= f; rg.exceptionvalue[3, s] *= f
    end
    return nothing
end
