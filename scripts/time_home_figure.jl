#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight

const OUT_PATH = joinpath(REPO_ROOT, "docs", "src", "assets", "coffee_scene_light_interception.png")
const BG = RGBf(0.982, 0.978, 0.969)
const HOME_FIGURE_PLOT_PAVING = 2500
const HOME_FIGURE_GROUND_SPAN = 45.0

sim_dir = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(sim_dir)
    cfg = read_light_config(sim_dir)
    scene = read_scene(cfg.scene)
    meteo = read_meteo(cfg.meteo)
    step = run_light_step(scene, first(meteo.rows), cfg)
    return scene, cfg, step
end

@benchmark run_archimed(sim_dir) # 25.237 s, with a memory estimate of 17.86 GiB, over 477428123 allocations.
