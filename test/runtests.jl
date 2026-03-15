using Test
using Artifacts
using ArchimedLight
using CairoMakie
using CSV
using Dates
using GeometryBasics
using LinearAlgebra: norm, cross
import PlantGeom
using ReferenceTests
using StaticArrays: SVector
using Tables

const _TEST_PROFILE = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_PROFILE", "all"))
const _RUN_RELEASE_TESTS = _TEST_PROFILE == "release"

@testset "Core tests" begin
    include("test_core.jl")
end

@testset "Fast fixtures" begin
    include("test_fast_fixtures.jl")
end

@testset "Synthetic scenes" begin
    include("synthetic_scene_cases.jl")
end

@testset "Model IO" begin
    include("test_model_io.jl")
end

if _RUN_RELEASE_TESTS
    @testset "Release tests" begin
        include("test_release.jl")
    end
end
