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
config_path = joinpath(REPO_ROOT, "example_2", "config.yml")
options, scene, meteo, models = read_config(config_path; plot_paving_override=HOME_FIGURE_PLOT_PAVING)

row = first(prepare_meteo(meteo, options).rows)
step = run_light_step(scene, models, row, options)
attach_light_step!(scene, step; fields=[:incident_par_flux])

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

mkpath(dirname(OUT_PATH))
save(OUT_PATH, fig)
println("wrote ", OUT_PATH)
