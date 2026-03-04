using ArchimedLight
using Documenter

DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true)

include("generate_reference_figures.jl")
generate_reference_figures()

makedocs(;
    modules=[ArchimedLight],
    authors="Rémi Vezy <VEZY@users.noreply.github.com> and contributors",
    sitename="ArchimedLight.jl",
    format=Documenter.HTML(;
        canonical="https://VEZY.github.io/ArchimedLight.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "ARCHIMED Reference" => "archimed_reference.md",
        "Full Example" => "full_example.md",
        "Composable Stages" => "stages.md",
        "API Reference" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/VEZY/ArchimedLight.jl",
    devbranch="main",
)
