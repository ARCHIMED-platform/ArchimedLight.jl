#!/usr/bin/env julia

using BenchmarkTools
const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight

function _load_config_bundle(config_path::AbstractString)
    raw = ArchimedLight._load_yaml_ordered(config_path)
    base = dirname(config_path)
    scene = read_scene(normpath(joinpath(base, string(raw["scene"]))))
    models = read_models(config_path)
    options = read_options(config_path)
    meteo = read_meteo(normpath(joinpath(base, string(raw["meteo"]))))
    return (scene=scene, models=models, options=options, meteo=meteo)
end

config_path = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(config_path)
    bundle = _load_config_bundle(config_path)
    row = first(prepare_meteo(bundle.meteo, bundle.options).rows)
    step = run_light_step(bundle.scene, bundle.models, row, bundle.options)
    return bundle.scene, bundle.models, bundle.options, step
end

@benchmark run_archimed(config_path)
# With java version:  java -jar ./example_1/archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar example_2/config.yml
# 2026-03-18 12:03:14.891 - Simulation time 14.448s0.0ms
# 2026-03-18 12:03:14.892 - maximum memory usage 2277,72 MB
