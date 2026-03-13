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

function _ri_par_colorrange(scene, step, cfg, window)
    meta = ArchimedLight._output_node_metadata(scene, cfg)
    vals = Float64[]
    ground_vals = Float64[]
    for nid in meta.node_ids
        v = get(step.budget.ri_par_f_per_node, nid, NaN)
        isfinite(v) || continue
        push!(vals, v)
        get(meta.group_per_node, nid, "") == "pavement" || continue
        x, y, _ = get(meta.barycenter_per_node, nid, (NaN, NaN, NaN))
        window.xmin <= x <= window.xmax || continue
        window.ymin <= y <= window.ymax || continue
        push!(ground_vals, v)
    end
    isempty(vals) && return (0.0, 1.0)
    if !isempty(ground_vals)
        hi = maximum(ground_vals)
        lo = max(0.0, hi - HOME_FIGURE_GROUND_SPAN)
        return (lo, hi)
    end
    hi = max(maximum(vals), eps(Float64))
    return (0.0, hi)
end


sim_dir = joinpath(REPO_ROOT, "example_2", "config.yml")
# cfg = read_light_config(_home_figure_config_path())
cfg = read_light_config(sim_dir)
scene = read_scene(cfg.scene)
meteo = read_meteo(cfg.meteo)
step = run_light_step(scene, first(meteo.rows), cfg)

colorrange = _ri_par_colorrange(scene, step, cfg, NamedTuple{(:xmin, :ymin, :xmax, :ymax)}(scene.scene_xy_bounds))
viz_mtg = visual_scene_mtg(
    scene,
    cfg,
    step;
    fields=[:ri_par_f_per_node],
    xy_bounds=scene.scene_xy_bounds,
)

CairoMakie.activate!()
fig, ax, p = plantviz(
    viz_mtg;
    color=:Ri_PAR_f,
    colormap=:viridis,
    colorrange=colorrange,
    figure=(size=(980, 700), backgroundcolor=BG),
)
ax.show_axis[] = false
PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")

mkpath(dirname(OUT_PATH))
save(OUT_PATH, fig)
println("wrote ", OUT_PATH)
