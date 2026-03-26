using AirspeedVelocity
using ArchimedLight
using BenchmarkTools
using CairoMakie
using Dates
using GeometryBasics
using StaticArrays: SVector
import LinearAlgebra: cross, norm

const PKG_ROOT = dirname(dirname(pathof(ArchimedLight)))
const FAST_FIXTURE_ROOT = joinpath(PKG_ROOT, "test", "fast_fixtures")
const SUITE = BenchmarkGroup()

include(joinpath(PKG_ROOT, "scripts", "generate_home_figure.jl"))

function _fixture_paths(name::String)
    root = joinpath(FAST_FIXTURE_ROOT, name, "input")
    return (
        root=root,
        config=joinpath(root, "config.yml"),
        meteo=joinpath(root, "meteo.csv"),
    )
end

function _load_fixture(name::String)
    paths = _fixture_paths(name)
    options, scene, meteo, models = ArchimedLight.read_config(paths.config)
    rows = ArchimedLight.prepare_meteo(meteo, options).rows
    return (paths=paths, scene=scene, models=models, options=options, meteo=meteo, rows=rows)
end

function _override_options(options0::ArchimedLight.LightOptions; kwargs...)
    params = Dict{Symbol,Any}()
    for (k, v) in kwargs
        if k == :pixel_size_m
            params[:pixel_size] = Float64(v)
        else
            params[k] = v
        end
    end
    return ArchimedLight.LightOptions(options0; params...)
end

function _synthetic_quad_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float64}[]
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    nodes = Dict{Int,ArchimedLight.SceneNodeData{Float64}}()
    xs = Float64[]
    ys = Float64[]

    for (i, spec) in enumerate(specs)
        p1 = ntuple(j -> Float64(spec.p1[j]), 3)
        p2 = ntuple(j -> Float64(spec.p2[j]), 3)
        p3 = ntuple(j -> Float64(spec.p3[j]), 3)
        p4 = ntuple(j -> Float64(spec.p4[j]), 3)
        append!(xs, (p1[1], p2[1], p3[1], p4[1]))
        append!(ys, (p1[2], p2[2], p3[2], p4[2]))

        base = length(points)
        append!(
            points,
            GeometryBasics.Point{3,Float64}[
                GeometryBasics.Point{3,Float64}(p1[1], p1[2], p1[3]),
                GeometryBasics.Point{3,Float64}(p2[1], p2[2], p2[3]),
                GeometryBasics.Point{3,Float64}(p3[1], p3[2], p3[3]),
                GeometryBasics.Point{3,Float64}(p4[1], p4[2], p4[3]),
            ],
        )
        append!(faces, GeometryBasics.TriangleFace{Int}[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, [i, i])

        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        area = area1 + area2
        barycenter = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            (p1[3] + p2[3] + p3[3] + p4[3]) / 4,
        )
        source_topology_id = Int(get(spec, :source_topology_id, i))
        object_id = Int(get(spec, :object_id, source_topology_id))
        group = String(get(spec, :group, "plate"))
        type_name = String(get(spec, :type, "plate"))
        nodes[i] = ArchimedLight.SceneNodeData(area, barycenter, group, type_name, source_topology_id, object_id)
    end

    return ArchimedLight.SceneGeometry(
        nothing,
        GeometryBasics.Mesh(points, faces),
        face2node,
        nodes,
        "benchmark_synthetic_scene",
        (minimum(xs), minimum(ys), maximum(xs), maximum(ys)),
    )
end

function _synthetic_meteo_row(;
    date::Dates.Date=Dates.Date(2020, 6, 21),
    start_time::Dates.Time=Dates.Time(12),
    duration_seconds::Float64=1.0,
    ri_par_f::Float64=100.0,
    ri_nir_f::Float64=0.0,
    direct_fraction::Float64=1.0,
    sun_azimut::Float64=180.0,
    sun_elevation::Float64=90.0,
)
    start_dt = Dates.DateTime(date, start_time)
    end_dt = start_dt + Dates.Millisecond(round(Int, duration_seconds * 1000))
    return (
        date=date,
        hour_start=Dates.format(Dates.Time(start_dt), Dates.DateFormat("HH:MM:SS")),
        hour_end=Dates.format(Dates.Time(end_dt), Dates.DateFormat("HH:MM:SS")),
        step_duration=duration_seconds,
        latitude=0.0,
        relativeHumidity=60.0,
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
        use="relativeHumidity RI_PAR_f",
    )
end

function _default_synthetic_models()
    ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=ArchimedLight.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
end

function _synthetic_fixture()
    scene = _synthetic_quad_scene([
        (
            p1=(0.0, 0.0, 1.0),
            p2=(1.0, 0.0, 1.0),
            p3=(1.0, 1.0, 1.0),
            p4=(0.0, 1.0, 1.0),
            group="upper",
            type="plate",
        ),
        (
            p1=(0.0, 0.0, 0.1),
            p2=(1.0, 0.0, 0.1),
            p3=(1.0, 1.0, 0.1),
            p4=(0.0, 1.0, 0.1),
            group="lower",
            type="plate",
        ),
    ])
    models = _default_synthetic_models()
    options = ArchimedLight.LightOptions()
    meteo = _synthetic_meteo_row(; ri_par_f=100.0, ri_nir_f=50.0)
    return (scene=scene, models=models, options=options, meteo=meteo)
end

const SIMPLEPLANT = _load_fixture("simpleplant_16_notoric")
const SKY_DIRECT = _load_fixture("sky_46_direct")
const SYNTHETIC = _synthetic_fixture()

SUITE["IO"] = BenchmarkGroup()
SUITE["IO"]["read config"]["simpleplant"] = @benchmarkable ArchimedLight.read_config($(SIMPLEPLANT.paths.config))
SUITE["IO"]["read models"]["simpleplant"] = @benchmarkable ArchimedLight.read_models($(SIMPLEPLANT.paths.config))
SUITE["IO"]["read options"]["simpleplant"] = @benchmarkable ArchimedLight.read_options($(SIMPLEPLANT.paths.config))
SUITE["IO"]["read scene"]["simpleplant"] = @benchmarkable ArchimedLight.read_scene(joinpath($(SIMPLEPLANT.paths.root), "scene", "simple.ops")) evals = 1
SUITE["IO"]["read meteo"]["simpleplant"] = @benchmarkable ArchimedLight.read_meteo($(SIMPLEPLANT.paths.meteo))

SUITE["Run light"] = BenchmarkGroup()
SUITE["Run light"]["step"] = BenchmarkGroup()
SUITE["Run light"]["series"] = BenchmarkGroup()

SUITE["Run light"]["step"]["simpleplant default"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, row, options) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        options = SIMPLEPLANT.options;
        row = first(SIMPLEPLANT.rows)
    ) evals = 1

SUITE["Run light"]["step"]["simpleplant toric + coarse pixels"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, row, options) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        options = _override_options(SIMPLEPLANT.options; toricity=true, pixel_size_m=0.40);
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, options).rows)
    ) evals = 1

SUITE["Run light"]["step"]["simpleplant scattering links"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, row, options; scattering_mode=:links) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        options = _override_options(SIMPLEPLANT.options; scattering=true, pixel_size_m=0.05);
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, options).rows)
    ) evals = 1

SUITE["Run light"]["step"]["sky direct fixture"] =
    @benchmarkable begin
        sky = ArchimedLight.compute_sky(row, options)
        turtle = ArchimedLight.build_turtle(options, sky)
        ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    end setup = (
        options = SKY_DIRECT.options;
        row = first(SKY_DIRECT.rows)
    ) evals = 1

SUITE["Run light"]["series"]["simpleplant uncached"] =
    @benchmarkable ArchimedLight.run_light_series(scene, models, meteo, options) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        meteo = SIMPLEPLANT.meteo;
        options = _override_options(SIMPLEPLANT.options; cache_radiation=false, scattering=false)
    ) evals = 1

SUITE["Run light"]["series"]["simpleplant cached"] =
    @benchmarkable ArchimedLight.run_light_series(scene, models, meteo, options) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        meteo = SIMPLEPLANT.meteo;
        options = _override_options(SIMPLEPLANT.options; cache_radiation=true, scattering=false)
    ) evals = 1

SUITE["Synthetic"] = BenchmarkGroup()
SUITE["Synthetic"]["first order stacked plates"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, meteo_row, options) setup = (
        scene = SYNTHETIC.scene;
        models = SYNTHETIC.models;
        meteo_row = SYNTHETIC.meteo;
        options = ArchimedLight.LightOptions(turtle_sectors=46, all_in_turtle=false, scattering=false, pixel_size=0.01)
    ) evals = 1

SUITE["Synthetic"]["stacked plates with scattering"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, meteo_row, options; scattering_mode=:raycast) setup = (
        scene = SYNTHETIC.scene;
        models = SYNTHETIC.models;
        meteo_row = SYNTHETIC.meteo;
        options = ArchimedLight.LightOptions(turtle_sectors=46, all_in_turtle=false, scattering=true, pixel_size=0.01)
    ) evals = 1

SUITE["Docs"] = BenchmarkGroup()
SUITE["Docs"]["home figure build"] =
    @benchmarkable simulate_home_figure(; repo_root=PKG_ROOT) evals = 1

function _leaf_benchmarks(group, prefix=String[])
    out = String[]
    for (name, item) in group
        next_prefix = [prefix; String(name)]
        if item isa BenchmarkGroup
            append!(out, _leaf_benchmarks(item, next_prefix))
        else
            push!(out, join(next_prefix, " / "))
        end
    end
    return out
end

if abspath(PROGRAM_FILE) == @__FILE__
    benches = sort(_leaf_benchmarks(SUITE))
    println("Loaded ArchimedLight benchmark suite for AirspeedVelocity/benchpkg.")
    println("Benchmark count: ", length(benches))
    for name in benches
        println(" - ", name)
    end
end
