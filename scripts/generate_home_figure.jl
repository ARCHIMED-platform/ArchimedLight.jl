#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight
using CairoMakie
using GeometryBasics
using PlantGeom

const OUT_PATH = joinpath(REPO_ROOT, "docs", "src", "assets", "coffee_scene_light_interception.png")
const BG = RGBf(0.982, 0.978, 0.969)
const HOME_FIGURE_PLOT_PAVING = 2500
const HOME_FIGURE_GROUND_SPAN = 45.0

function _coffee_fixture()
    cfg = read_light_config(
        joinpath(REPO_ROOT, "java_implementation", "archimed-lib-2018", "tests", "test-cafeier", "config.yml"),
    )
    cfg.raw["plot_paving_override"] = HOME_FIGURE_PLOT_PAVING
    scene = read_scene(cfg.scene)
    meteo = read_meteo(cfg.meteo)
    step = run_light_step(scene, first(meteo.rows), cfg)
    return cfg, scene, step
end

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

function _local_window(scene; pad_factor=0.18, min_pad=0.45)
    plant_xyz = collect(values(scene.barycenter_per_node))
    xs = [xyz[1] for xyz in plant_xyz if isfinite(xyz[1])]
    ys = [xyz[2] for xyz in plant_xyz if isfinite(xyz[2])]
    zs = [xyz[3] for xyz in plant_xyz if isfinite(xyz[3])]
    isempty(xs) && error("Scene has no finite node barycenters.")
    isempty(ys) && error("Scene has no finite node barycenters.")
    isempty(zs) && error("Scene has no finite node barycenters.")
    xpad = max(min_pad, pad_factor * (maximum(xs) - minimum(xs)))
    ypad = max(min_pad, pad_factor * (maximum(ys) - minimum(ys)))
    zmin = minimum(zs)
    zmax = maximum(zs)
    xmin = minimum(xs) - xpad
    xmax = maximum(xs) + xpad
    ymin = minimum(ys) - ypad
    ymax = maximum(ys) + ypad
    return (; xmin, xmax, ymin, ymax, zmin, zmax)
end

function _display_window(scene)
    plant = _local_window(scene)
    xmax = max(plant.xmax, plant.xmin + 4.8)
    ymax = max(plant.ymax, plant.ymin + 5.8)
    return (; xmin=plant.xmin, xmax, ymin=plant.ymin, ymax, zmin=plant.zmin, zmax=plant.zmax)
end

function _set_camera!(ax, scene)
    window = _display_window(scene)
    xmid = (window.xmin + window.xmax) / 2
    ymid = (window.ymin + window.ymax) / 2
    zmid = 0.22 * window.zmax
    xy_span = max(window.xmax - window.xmin, window.ymax - window.ymin)
    z_span = max(window.zmax - window.zmin, 0.6)
    radius = max(1.2, 0.55 * max(xy_span, z_span))
    eye = Point3f(xmid + 1.18radius, ymid - 0.78radius, zmid + 0.96radius)
    lookat = Point3f(xmid, ymid, zmid)
    update_cam!(ax.scene, eye, lookat, Vec3f(0, 0, 1))
    ax.scene.camera_controls.fov[] = 19.0
    return nothing
end

function main()
    cfg, scene, step = _coffee_fixture()
    window = _display_window(scene)
    colorrange = _ri_par_colorrange(scene, step, cfg, window)
    viz_mtg = visual_scene_mtg(
        scene,
        cfg,
        step;
        fields=[:ri_par_f_per_node],
        xy_bounds=(window.xmin, window.ymin, window.xmax, window.ymax),
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
    _set_camera!(ax, scene)
    PlantGeom.colorbar(fig[1, 2], p, label="Ri_PAR_f (W m^-2)")

    mkpath(dirname(OUT_PATH))
    save(OUT_PATH, fig)
    println("wrote ", OUT_PATH)
end

main()
