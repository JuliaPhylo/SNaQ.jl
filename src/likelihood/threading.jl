# Per-thread scratch space for the loss and gradient.

"""
    useparallelloop(n)

Whether a loop over `n` independent quartets is worth spreading across the thread pool.
Dispatching to worker threads costs the same whether or not there is much work to spread,
and these loops run once per gradient evaluation, so small `n` is faster served serially.
"""
@inline useparallelloop(n::Int)::Bool = Threads.nthreads() > 1 && n >= 4 * Threads.nthreads()


# Scratch buffers, keyed by parameter count, with one column (or one entry) per thread.
const THREAD_QUARTET_GRAD::Dict{Int, Array{Float64, 3}} = Dict{Int, Array{Float64, 3}}()
const THREAD_TOTAL_GRAD::Dict{Int, Matrix{Float64}} = Dict{Int, Matrix{Float64}}()

# Each thread needs its OWN `RunningGradient`: `params_seen` must not be packed into bits
# shared with another thread's, or the two read-modify-write the same 64-bit word and
# silently clobber each other's flags.
const THREAD_RUNNING_GRAD::Dict{Int, Vector{RunningGradient}} = Dict{Int, Vector{RunningGradient}}()

# Thread id => column index, skipping `:foreign` threads. Filled by `initthreadscratch!`
# rather than at compile time: package images are built single-threaded, so a baked-in
# value has length 1 and every thread above the first indexes out of bounds under `-t N`.
const THREAD_COLUMN_REF::Base.RefValue{Vector{Int}} = Ref(Int[])

@inline threadcolumns()::Vector{Int} = THREAD_COLUMN_REF[]


"""
    initthreadscratch!()

Sizes the thread-indexed scratch to the current process's thread count, discarding anything
carried over from the session that built the package image. Called from `SNaQ.__init__`.
"""
function initthreadscratch!()
    m = zeros(Int, Threads.maxthreadid())
    j = 1
    for tid = 1:Threads.maxthreadid()
        if Threads.threadpool(tid) != :foreign
            m[tid] = j
            j += 1
        end
    end
    THREAD_COLUMN_REF[] = m
    empty!(THREAD_QUARTET_GRAD)
    empty!(THREAD_RUNNING_GRAD)
    empty!(THREAD_TOTAL_GRAD)
    return nothing
end


"""
    threadscratch(nparam)

The per-quartet gradient buffer, [`RunningGradient`](@ref)s and total-gradient buffer for
`nparam` parameters, created on first use. Returned as a tuple, one column per thread.
"""
function threadscratch(nparam::Int)
    nthr = Threads.maxthreadid()
    # `Threads.maxthreadid()` can grow after `__init__` (a C library adopting a `:foreign`
    # thread, say); resize here, where no parallel region is running.
    length(threadcolumns()) < nthr && initthreadscratch!()
    if !haskey(THREAD_TOTAL_GRAD, nparam)
        THREAD_QUARTET_GRAD[nparam] = zeros(nparam, 3, nthr)
        THREAD_TOTAL_GRAD[nparam] = zeros(nparam, nthr)
        THREAD_RUNNING_GRAD[nparam] = [RunningGradient(nparam) for _ = 1:nthr]
    end
    return THREAD_QUARTET_GRAD[nparam], THREAD_RUNNING_GRAD[nparam], THREAD_TOTAL_GRAD[nparam]
end
