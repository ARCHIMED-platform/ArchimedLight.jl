#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
if abspath(PROGRAM_FILE) == @__FILE__
    import Pkg
    Pkg.activate(joinpath(REPO_ROOT, "test"))
end

using ArchimedLight
using CairoMakie
using PlantGeom

const OUT_PATH = joinpath(REPO_ROOT, "docs", "src", "assets", "coffee_scene_light_interception.png")
const BG = RGBf(0.982, 0.978, 0.969)
const HOME_FIGURE_PLOT_PAVING = 2500

function simulate_home_figure(
    repo_root::AbstractString=REPO_ROOT,
    plot_paving_override::Int=HOME_FIGURE_PLOT_PAVING,
)
    config_path = joinpath(repo_root, "example_2", "config.yml")
    options, scene, meteo, models = ArchimedLight.read_config(config_path; plot_paving_override=plot_paving_override)
    row = first(ArchimedLight.prepare_meteo(meteo, options).rows)
    step = ArchimedLight.run_light_step(scene, models, row, options)
    return scene, models, options, step
end

function build_home_figure(;
    repo_root::AbstractString=REPO_ROOT,
    plot_paving_override::Int=HOME_FIGURE_PLOT_PAVING,
)
    scene, _, _, step = simulate_home_figure(repo_root, plot_paving_override)
    ArchimedLight.attach_light_step!(scene, step; fields=[:incident_par_flux])
    CairoMakie.activate!()
    fig, ax, p = plantviz(
        scene.mtg;
        color=:Ri_PAR_f,
        colormap=:thermal,
        colorrange=(0.0, 450.0),
        figure=(size=(980, 700), backgroundcolor=BG),
    )
    ax.show_axis[] = false
    PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")
    return fig
end

function generate_home_figure(out_path::AbstractString=OUT_PATH; kwargs...)
    fig = build_home_figure(; kwargs...)
    mkpath(dirname(out_path))
    save(out_path, fig)
    println("wrote ", out_path)
    return out_path
end

if abspath(PROGRAM_FILE) == @__FILE__
    generate_home_figure()
end
