#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight
using CairoMakie
using GeometryBasics
using OrderedCollections: OrderedDict
using PlantGeom

const OUT_PATH = joinpath(REPO_ROOT, "docs", "src", "assets", "coffee_scene_light_interception.png")
const BG = RGBf(0.982, 0.978, 0.969)
const HOME_FIGURE_PLOT_PAVING = 2500
const HOME_FIGURE_GROUND_SPAN = 45.0

sim_dir = joinpath(REPO_ROOT, "example_2", "config.yml")
cfg = read_light_config(sim_dir)
cfg.models["pavement"].types["Cobblestone"].plot_paving = HOME_FIGURE_PLOT_PAVING
scene = read_scene(cfg.paths.scene)
meteo = read_meteo(cfg.paths.meteo)
step = run_light_step(scene, first(meteo.rows), cfg)

viz_mtg = visual_scene_mtg(
    scene,
    cfg,
    step;
    fields=[:incident_par_flux],
    xy_bounds=scene.scene_xy_bounds,
)

CairoMakie.activate!()
fig, ax, p = plantviz(
    viz_mtg;
    color=:Ri_PAR_f,
    colormap=:thermal,#, :batlow, :thermal, :cividis
    # colorrange=colorrange,
    colorrange=(0.0, 450.0),
    figure=(size=(980, 700), backgroundcolor=BG),
)
ax.show_axis[] = false
PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")

mkpath(dirname(OUT_PATH))
save(OUT_PATH, fig)
println("wrote ", OUT_PATH)
