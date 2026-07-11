#!/usr/bin/env julia

using BenchmarkTools
const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight

config_path = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(config_path)
    options, scene, meteo, models = read_config(config_path)
    row = first(prepare_meteo(meteo, options))
    step = run_light_step(scene, models, row, options)
    return scene, models, options, step
end

trial = @benchmark run_archimed($config_path)
display(trial)

# With Julia: 
# Single result which took 1.6243 s (12.79%) to evaluate,
#  with a memory estimate of 1.70 GiB, allocs estimate: 12139919.
