#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight

sim_dir = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(sim_dir)
    cfg = read_light_config(sim_dir)
    scene = read_scene(cfg.paths.scene)
    meteo = read_meteo(cfg.paths.meteo)
    step = run_light_step(scene, first(meteo.rows), cfg)
    return scene, cfg, step
end

@benchmark run_archimed(sim_dir) # 25.237 s, with a memory estimate of 17.86 GiB, over 477428123 allocations.
# With java version:  java -jar ./example_1/archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar example_2/config.yml (but first, change meteo step into 1,1 in the config file)
# 2026-03-18 12:03:14.891 - Simulation time 14.448s0.0ms
# 2026-03-18 12:03:14.892 - maximum memory usage 2277,72 MB