using SNaQ
using Documenter

DocMeta.setdocmeta!(SNaQ, :DocTestSetup, :(using SNaQ); recursive=true)

makedocs(;
    modules=[SNaQ],
    authors="Claudia Solis-Lemus <crsl4@users.noreply.github.com>, Paul Bastide <pbastide@users.noreply.github.com>, Cecile Ane <cecileane@users.noreply.github.com>, and contributors",
    sitename="SNaQ.jl",
    format=Documenter.HTML(;
        canonical="https://JuliaPhylo.github.io/SNaQ.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaPhylo/SNaQ.jl",
    devbranch="master",
)
