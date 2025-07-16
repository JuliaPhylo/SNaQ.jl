using Pkg
Pkg.develop(PackageSpec(path = joinpath(@__DIR__, "..")))
Pkg.instantiate()


using SNaQ, PhyloNetworks
using Documenter


##add note to all pages in documentation
note = """
    !!! info "Important Note:"
        This documentation pertains to SNaQ v1.1 and may differ from the specific implementation
        originally described in [Solís-Lemus & Ané (2016)](https://doi.org/10.1371/journal.pgen.1005896).
        See documentation SNaQ v1.0  for the original implementation.
    """

const original_src = joinpath(@__DIR__, "src")
const temp_src = joinpath(@__DIR__, "temp_src_for_build") # temp directory with note added to md files

# Ensure the temporary directory is clean and exists
if isdir(temp_src)
    rm(temp_src; recursive=true, force=true)
end
mkpath(temp_src)
for (root, dirs, files) in walkdir(original_src) #go thru src directory
    relative_path = relpath(root, original_src)  #Calculate the relative path from the original src directory
    
    temp_root_path = joinpath(temp_src, relative_path) # Determine the corresponding path in the temporary directory
    
    # Create directories in the temporary structure
    mkpath(temp_root_path)

    for file in files
        original_filepath = joinpath(root, file)
        temp_filepath = joinpath(temp_root_path, file)

        if endswith(file, ".md") ## add note to md files
            content = read(original_filepath, String)
            # Prepend the note
            new_content = note * "\n" * content
            write(temp_filepath, new_content)
        else #copy file directly
            cp(original_filepath, temp_filepath; force=true)
        end
    end
end


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
            "Improving runtimes" => "man/parallelcomputation.md",
            "Multiple alleles" => "man/multiplealleles.md",
            "Error reporting" => "man/error_reporting.md"
        ],
        "Library" => [
            "Public" => "lib/public.md",
            "Internals" => "lib/internals.md",
        ]
    ],
    source= temp_src,
    plugins = [links]
)

deploydocs(;
    repo="github.com/JuliaPhylo/SNaQ.jl",
    push_preview = true,
    devbranch="main",
)

# Clean up the temporary directory
rm(temp_src; recursive=true, force=true)