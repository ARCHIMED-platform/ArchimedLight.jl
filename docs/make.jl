using ArchimedLight
using Documenter

# For interactive rendering:
# using LiveServer
# servedocs()
# or julia --project=docs -e 'using LiveServer; LiveServer.serve(dir = "docs/build/1")'

DocMeta.setdocmeta!(ArchimedLight, :DocTestSetup, :(using ArchimedLight); recursive=true)

# Update the videos:
# include(joinpath(@__DIR__, "make_video", "make_video.jl"))

makedocs(;
    modules=[ArchimedLight],
    authors="Rémi Vezy <VEZY@users.noreply.github.com> and contributors",
    repo=Documenter.Remotes.GitHub("VEZY", "ArchimedLight.jl"),
    sitename="ArchimedLight.jl",
    format=Documenter.HTML(
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://VEZY.github.io/ArchimedLight.jl",
        edit_link="main",
        assets=String[],
        size_threshold=700000,
    ),
    pages=[
        "Home" => "index.md",
        "Getting Started" => "getting_started.md",
        "Tutorials" => [
            "Beginner Workflows" => "tutorial_beginner_workflows.md",
            "File-Based Workflow" => "tutorial_files.md",
            "Interactive Workflow" => "tutorial_interactive.md",
            "Copy-Paste Templates" => "onboarding_templates.md",
        ],
        "Reference" => [
            "Configuration And Options" => "reference_config.md",
            "Scene Files And Semantics" => "reference_scene.md",
            "Model Files" => "reference_models.md",
            "Meteo Inputs" => "reference_meteo.md",
            "Voxel Light Interception" => "voxel_interception.md",
            "Voxel Terrain And Soil" => "voxel_terrain.md",
        ],
        "Outputs" => "outputs.md",
        "Performance Benchmarks" => "performance_benchmarks.md",
        "Under The Hood" => [
            "Pipeline Overview" => "theory_pipeline.md",
            "First-Order Interception" => "theory_interception.md",
            "Scattering And Assumptions" => "theory_scattering.md",
            "Voxel Scattering Theory" => "theory_voxel_scattering.md",
            "Voxel Future Directions" => "voxel_future_directions.md",
        ],
        "Advanced Usage" => [
            "File-Based Example" => "full_example.md",
            "Composable Stages" => "stages.md",
            "Historical ARCHIMED Reference" => "archimed_reference.md",
        ],
        "API Reference" => [
            "Public API" => "api.md",
            "Advanced API" => "advanced_api.md",
        ],
    ],
)

include(joinpath(@__DIR__, "check_inline_voxel_plot.jl"))
check_inline_voxel_plot(joinpath(@__DIR__, "build"))

deploydocs(; repo="github.com/VEZY/ArchimedLight.jl.git", devbranch="main", push_preview=true)
# Visit https://VEZY.github.io/ArchimedLight.jl/previews/PR26 to visualize the preview of the PR #26
