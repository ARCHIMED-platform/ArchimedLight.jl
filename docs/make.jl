using ArchimedLight
using Documenter

DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/VEZY/ArchimedLight.jl",
    devbranch="main",
)
