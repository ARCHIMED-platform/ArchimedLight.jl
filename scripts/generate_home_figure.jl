#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight
using CairoMakie
using PlantGeom

const OUT_PATH = joinpath(REPO_ROOT, "docs", "src", "assets", "coffee_scene_light_interception.png")
const BG = RGBf(0.982, 0.978, 0.969)
const HOME_FIGURE_PLOT_PAVING = 2500

function _load_config_bundle(config_path::AbstractString)
    raw = ArchimedLight._load_yaml_ordered(config_path)
    base = dirname(config_path)
    scene = read_scene(normpath(joinpath(base, string(raw["scene"]))))
    models = read_models(config_path)
    options = read_options(config_path)
    meteo = read_meteo(normpath(joinpath(base, string(raw["meteo"]))))
    return (scene=scene, models=models, options=options, meteo=meteo)
end

bundle = _load_config_bundle(joinpath(REPO_ROOT, "example_2", "config.yml"))
side = round(Int, sqrt(HOME_FIGURE_PLOT_PAVING))
add_ground!(bundle.scene; nx=side, ny=side, group="pavement", type="Cobblestone")

row = first(prepare_meteo(bundle.meteo, bundle.options).rows)
step = run_light_step(bundle.scene, bundle.models, row, bundle.options)
attach_light_step!(bundle.scene, step; fields=[:incident_par_flux])

CairoMakie.activate!()
fig, ax, p = plantviz(
    bundle.scene.mtg;
    color=:Ri_PAR_f,
    colormap=:thermal,
    colorrange=(0.0, 450.0),
    figure=(size=(980, 700), backgroundcolor=BG),
)
ax.show_axis[] = false
PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")

mkpath(dirname(OUT_PATH))
save(OUT_PATH, fig)
println("wrote ", OUT_PATH)
