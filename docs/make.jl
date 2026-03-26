using ArchimedLight
using Documenter

DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive = true)

include("generate_reference_figures.jl")
generate_reference_figures()

makedocs(;
    modules = [ArchimedLight],
    authors = "Rémi Vezy <VEZY@users.noreply.github.com> and contributors",
    repo = "https://github.com/VEZY/ArchimedLight.jl/blob/{commit}{path}#{line}",
    sitename = "ArchimedLight.jl",
    format = Documenter.HTML(; canonical = "https://VEZY.github.io/ArchimedLight.jl"),
    pages = [
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => [
            "File-Based Workflow" => "tutorial_files.md",
            "Interactive Workflow" => "tutorial_interactive.md",
        ],
        "Reference" => [
            "Configuration And Options" => "reference_config.md",
            "Scene Files And Semantics" => "reference_scene.md",
            "Model Files" => "reference_models.md",
            "Meteo Inputs" => "reference_meteo.md",
        ],
        "Outputs" => "outputs.md",
        "Under The Hood" => [
            "Pipeline Overview" => "theory_pipeline.md",
            "First-Order Interception" => "theory_interception.md",
            "Scattering And Assumptions" => "theory_scattering.md",
        ],
        "Advanced Usage" => [
            "Full Example" => "full_example.md",
            "Composable Stages" => "stages.md",
            "Historical ARCHIMED Reference" => "archimed_reference.md",
        ],
        "API Reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/VEZY/ArchimedLight.jl")
