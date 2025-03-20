using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()


using SNaQ, PhyloNetworks
using Documenter

# Interlink with PhyloNetworks
using DocumenterInterLinks
links = InterLinks(
    "PhyloNetworks" => "https://juliaphylo.github.io/PhyloNetworks.jl/stable/"
)

# NOTE: default loading of PhyloNetworks in all docstring examples
DocMeta.setdocmeta!(SNaQ, :DocTestSetup, :(using PhyloNetworks, SNaQ); recursive=true)

makedocs(;
    modules=[SNaQ],
    authors="Claudia Solis-Lemus <crsl4@users.noreply.github.com>, Cécile Ané <cecileane@users.noreply.github.com>, and contributors",
    sitename="SNaQ.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaPhylo.github.io/SNaQ.jl",
        edit_link="dev",
        assets=String[],
        prettyurls = get(ENV, "CI", nothing) == "true", # easier local build
        size_threshold = 600 * 2^10,
        size_threshold_warn = 500 * 2^10, # 600 KiB
    ),
    # exception, so warning-only for :missing_docs. List all others:
    warnonly = Documenter.except(:autodocs_block, :cross_references, :docs_block,
        :doctest, :eval_block, :example_block, :footnote, :linkcheck_remotes,
        :linkcheck, :meta_block, :parse_error, :setup_block),
    pages=[
        "Home" => "index.md",
        "Manual" => [
            "Installation" => "man/installation.md",
            "Network estimation" => "man/snaq_est.md",
            "Candidate networks" => "man/fixednetworkoptim.md",
            "Extract expected CFs" => "man/expectedCFs.md",
            "Bootstrap" => "man/bootstrap.md",
            "Parallel computation" => "man/parallelcomputation.md",
            "Multiple alleles" => "man/multiplealleles.md",
            "Error reporting" => "man/error_reporting.md"
        ],
        "Library" => [
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
    plugins = [links]
)

deploydocs(;
    repo="github.com/JuliaPhylo/SNaQ.jl",
    push_preview = true,
    devbranch="main",
)
