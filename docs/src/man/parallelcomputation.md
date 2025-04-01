# Parallel computations

> This documentation pertains to SNaQ v1.0 as originally described in [Solís-Lemus & Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)

For network estimation, multiple runs can done in parallel.
For example, if your machine has 4 or more processors (or cores),
you can tell julia to use 4 processors by starting julia with `julia -p 4`,
or by starting julia the usual way (`julia`) and then adding processors with:

```julia
using Distributed
addprocs(4)
```

If we load a package (`using SNaQ`) before adding processors,
then we need to re-load it again so that all processors have access to it:

```julia
@everywhere using PhyloNetworks, SNaQ
```

After that, running any of the `snaq!(...)` command will use
different cores for different runs, as processors become available.
Fewer details are printed to the log file when multiple cores
are used in parallel.

When running [`bootsnaq`](@ref), the analysis of each bootstrap replicate
will use multiple cores to parallelize separate runs of that particular
bootstrap replicate. You may parallelize things further by running
`bootsnaq` multiple times (on separate machines for instance), each time
for a small subset of bootstrap replicates, and with a different seed each time.

We may tell julia to add more processors than our machine has,
but we will not receive any performance benefits.
At any time during the julia session, `nworkers()` tells us how many
worker processors julia has access to.

Below is an example of how to use a cluster, to run many independent
[`snaq!`](@ref) searches in parallel on a cluster running the
[slurm](https://slurm.schedmd.com) job manager
(other managers would require a different, but similar submit file).
This example uses 2 files:
1. a julia script file, to do many runs of `snaq!` in parallel,
   asking for many cores (default: 10 runs, asking for 10 cores).
   This julia script can take arguments: the maximum allowed
   number of hybridizations `hmax`, and the number of runs
   (to run 50 runs instead of 10, say).
2. a submit file, to launch the julia script.

**First**: the example julia script, below or [here](https://github.com/juliaphylo/SNaQ/blob/main/examples/runSNaQ.jl), is assumed (by the submit file)
to be called `runSNaQ.jl`. It uses a starting tree that
is assumed to be available in a file named `astraltree.tre`, but that
could be modified
(to use a network with h=1 to start the search with hmax=2 for instance).
It also assumes that the quartet concordance factor data are in file
`tableCF_speciesNames.csv`. Again, this file name should be adjusted.
To run this julia script for 50 runs and hmax=3, do `julia runSNaQ.jl 3 50`.

```julia
#!/usr/bin/env julia

# file "runSNaQ.jl". run in the shell like this in general:
# julia runSNaQ.jl hvalue nruns
# example for h=2 and default 10 runs:
# julia runSNaQ.jl 2
# or example for h=3 and 50 runs:
# julia runSNaQ.jl 3 50

length(ARGS) > 0 ||
    error("need 1 or 2 arguments: # reticulations (h) and # runs (optional, 10 by default)")
h = parse(Int, ARGS[1])
nruns = 10
if length(ARGS) > 1
    nruns = parse(Int, ARGS[2])
end
outputfile = string("net", h, "_", nruns, "runs") # example: "net2_10runs"
seed = 1234 + h # change as desired! Best to have it different for different h
@info "will run SNaQ with h=$h, # of runs=$nruns, seed=$seed, output will go to: $outputfile"

using Distributed
addprocs(nruns)
@everywhere using SNaQ
net0 = readnewick("astraltree.tre");
using DataFrames, CSV
df_sp = CSV.read("tableCF_speciesNames.csv", DataFrame; pool=false);
d_sp = readtableCF!(df_sp);
net = snaq!(net0, d_sp, hmax=h, filename=outputfile, seed=seed, runs=nruns)
```

When julia is called on a script, whatever comes after "julia scriptname"
is given to julia in an array of values. This array is called `ARGS`.
So if we call a script like this: `julia runSNaQ.jl 2`
then the script will know the arguments through `ARGS`,
which would contain a single element, `"2"`.
This first element is just a string, at this stage. We want to use it as a number,
so we need to ask julia to parse the string into an integer.

**Second**: we need a "submit" file to ask a job scheduler like
[slurm](https://slurm.schedmd.com) to submit our julia script to a cluster.
In the submit file below, the first 5 lines set things up for slurm.
They are most likely to be specific to your cluster.
The main idea here is to use a slurm "array" from 0 to 3, to run our
julia script multiple times, 4 times actually: from hmax=0 to hmax=3.
Each would do 30 runs
(and each would be allocated 30 cores in the submit script below).
Then log out of the cluster and go for coffee.

```bash
#!/bin/bash
#SBATCH -o path/to/slurm/log/file/runsnaq_slurm%a.log
#SBATCH -J runsnaq
#SBATCH --array=0-3
#SBATCH -c 30
## --array: to run multiple instances of this script,
##          one for each value in the array.
##          1 instance = 1 task
## -J job name
## -c number of cores (CPUs) per task

echo "slurm task ID = $SLURM_ARRAY_TASK_ID used as hmax"
echo "start of SNaQ parallel runs on $(hostname)"
# finally: launch the julia script, using Julia executable appropriate for slurm, with full paths:
/workspace/software/bin/julia --history-file=no -- runSNaQ.jl $SLURM_ARRAY_TASK_ID 30 > net${SLURM_ARRAY_TASK_ID}_30runs.screenlog 2>&1
echo "end of SNaQ run ..."
```

> This documentation pertains to SNaQ v1.0 as originally described in [Solís-Lemus & Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896)
