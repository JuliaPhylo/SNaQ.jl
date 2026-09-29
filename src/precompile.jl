# Module initialization and the precompilation workload.

using PrecompileTools

"""
Thread-indexed scratch must be sized in the process that runs, never baked into the package
image: images are built single-threaded (see `likelihood/threading.jl`).
"""
function __init__()
    initthreadscratch!()
    return nothing
end


# Bake the lazy `snaq!` call tree into the package image: without this the first call in a
# session spends ~9 seconds compiling before doing any work.
@setup_workload begin
    pc_start = readnewick("(((t1:1.0,t2:1.0):0.5,(t3:1.0,t4:1.0):0.5):0.5,(t5:1.0,t6:1.0):0.5);")
    pc_gts = [
        readnewick("(((t1:1.0,t2:1.0):0.5,(t3:1.0,t4:1.0):0.5):0.5,(t5:1.0,t6:1.0):0.5);"),
        readnewick("(((t1:1.0,t3:1.0):0.5,(t2:1.0,t4:1.0):0.5):0.5,(t5:1.0,t6:1.0):0.5);"),
        readnewick("(((t1:1.0,t2:1.0):0.5,(t3:1.0,t5:1.0):0.5):0.5,(t4:1.0,t6:1.0):0.5);"),
    ]
    @compile_workload begin
        redirect_stdout(devnull) do
            snaq!(pc_start, DataCF(pc_gts; lazy=true); hmax=0, propQuartets=0.5,
                  propQuartetsFinal=0.9, runs=1, Nfail=1, filename="", seed=1, verbose=false)
        end
    end
end
