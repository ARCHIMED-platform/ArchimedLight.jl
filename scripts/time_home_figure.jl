#!/usr/bin/env julia

using BenchmarkTools
const REPO_ROOT = dirname(@__DIR__)
import Pkg
Pkg.activate(joinpath(REPO_ROOT, "test"))

using ArchimedLight
# using KernelAbstractions
# using Metal
# metal_backend = KernelAbstractions.get_backend(MtlArray(zeros(Float32, 1)))
# interception = RaycoreInterceptionBackend(backend=metal_backend)

config_path = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(config_path)
    sim, meteo = read_simulation(
        config_path;
        # interception_backend=interception,
        # scattering_backend=RaycoreScatteringBackend(interception),
    )
    step = run_light(sim, first(meteo))
    return step
end

trial = @benchmark run_archimed($config_path)
display(trial)
# With java version:  java -jar ./example_1/archimed-lib-2018-0.0.1-SNAPSHOT-jar-with-dependencies.jar example_2/config.yml
# 2026-03-18 12:03:14.891 - Simulation time 14.448s0.0ms
# 2026-03-18 12:03:14.892 - maximum memory usage 2277,72 MB

# With Julia: 
# Single result which took 1.6243 s (12.79%) to evaluate,
#  with a memory estimate of 1.70 GiB, allocs estimate: 12139919.
