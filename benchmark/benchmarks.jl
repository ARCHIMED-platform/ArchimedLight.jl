using AirspeedVelocity
using ArchimedLight
using BenchmarkTools
using Dates
using GeometryBasics
using StaticArrays: SVector
import PlantGeom
import LinearAlgebra: cross, norm

const PKG_ROOT = dirname(dirname(pathof(ArchimedLight)))
const FAST_FIXTURE_ROOT = joinpath(PKG_ROOT, "test", "fast_fixtures")
const SUITE = BenchmarkGroup()

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
    cfg = ArchimedLight.read_light_config(paths.config)
    scene = ArchimedLight.read_scene(cfg.scene)
    meteo = ArchimedLight.read_meteo(cfg.meteo)
    rows = ArchimedLight.prepare_meteo(meteo, cfg).rows
    return (paths=paths, cfg=cfg, scene=scene, meteo=meteo, rows=rows)
end

function _override_cfg(cfg0::ArchimedLight.LightConfig; kwargs...)
    cfg = deepcopy(cfg0)
    for (k, v) in kwargs
        if k == :pixel_size_m
            cfg.raw["pixel_size"] = Float64(v) * 100.0
        else
            cfg.raw[String(k)] = v
        end
    end
    ArchimedLight.refresh_light_config!(cfg; reload_models=true)
    cfg.raw["output_directory"] = mktempdir()
    return cfg
end

function _synthetic_quad_scene(specs::AbstractVector{<:NamedTuple})
    points = GeometryBasics.Point{3,Float32}[]
    faces = PlantGeom.Face3[]
    face2node = Int[]
    total_area_per_node = Dict{Int,Float64}()
    barycenter_per_node = Dict{Int,NTuple{3,Float64}}()
    node_group = Dict{Int,String}()
    node_type = Dict{Int,String}()
    java_item_id_per_node = Dict{Int,Int}()
    java_component_id_per_node = Dict{Int,Int}()
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
            GeometryBasics.Point{3,Float32}[
                GeometryBasics.Point{3,Float32}(Float32(p1[1]), Float32(p1[2]), Float32(p1[3])),
                GeometryBasics.Point{3,Float32}(Float32(p2[1]), Float32(p2[2]), Float32(p2[3])),
                GeometryBasics.Point{3,Float32}(Float32(p3[1]), Float32(p3[2]), Float32(p3[3])),
                GeometryBasics.Point{3,Float32}(Float32(p4[1]), Float32(p4[2]), Float32(p4[3])),
            ],
        )
        append!(faces, PlantGeom.Face3[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, [i, i])

        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        total_area_per_node[i] = area1 + area2
        barycenter_per_node[i] = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            (p1[3] + p2[3] + p3[3] + p4[3]) / 4,
        )
        node_group[i] = String(get(spec, :group, "plate"))
        node_type[i] = String(get(spec, :type, "plate"))
        java_item_id_per_node[i] = Int(get(spec, :item_id, i))
        java_component_id_per_node[i] = Int(get(spec, :component_id, 1))
    end

    return ArchimedLight.SceneGeometry(
        nothing,
        GeometryBasics.Mesh(points, faces),
        face2node,
        total_area_per_node,
        barycenter_per_node,
        node_group,
        node_type,
        java_item_id_per_node,
        java_component_id_per_node,
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

function _synthetic_fixture()
    cfg_ref = ArchimedLight.read_light_config(joinpath(FAST_FIXTURE_ROOT, "simpleplant_16_notoric", "input", "config.yml"))
    scene = _synthetic_quad_scene([
        (
            p1=(0.0, 0.0, 1.0),
            p2=(1.0, 0.0, 1.0),
            p3=(1.0, 1.0, 1.0),
            p4=(0.0, 1.0, 1.0),
            item_id=1,
            component_id=1,
            group="upper",
            type="plate",
        ),
        (
            p1=(0.0, 0.0, 0.1),
            p2=(1.0, 0.0, 0.1),
            p3=(1.0, 1.0, 0.1),
            p4=(0.0, 1.0, 0.1),
            item_id=2,
            component_id=1,
            group="lower",
            type="plate",
        ),
    ])
    meteo = _synthetic_meteo_row(; ri_par_f=100.0, ri_nir_f=50.0)
    return (cfg_ref=cfg_ref, scene=scene, meteo=meteo)
end

const SIMPLEPLANT = _load_fixture("simpleplant_16_notoric")
const SKY_DIRECT = _load_fixture("sky_46_direct")
const SYNTHETIC = _synthetic_fixture()

SUITE["IO"] = BenchmarkGroup()
SUITE["IO"]["read config"]["simpleplant"] = @benchmarkable ArchimedLight.read_light_config($(SIMPLEPLANT.paths.config))
SUITE["IO"]["read scene"]["simpleplant"] = @benchmarkable ArchimedLight.read_scene($(SIMPLEPLANT.cfg.scene)) evals = 1
SUITE["IO"]["read meteo"]["simpleplant"] = @benchmarkable ArchimedLight.read_meteo($(SIMPLEPLANT.cfg.meteo))

SUITE["Run light"] = BenchmarkGroup()
SUITE["Run light"]["step"] = BenchmarkGroup()
SUITE["Run light"]["series"] = BenchmarkGroup()

SUITE["Run light"]["step"]["simpleplant default"] =
    @benchmarkable ArchimedLight.run_light_step(scene, row, cfg) setup = (
        scene = SIMPLEPLANT.scene;
        cfg = SIMPLEPLANT.cfg;
        row = first(SIMPLEPLANT.rows)
    ) evals = 1

SUITE["Run light"]["step"]["simpleplant toric + coarse pixels"] =
    @benchmarkable ArchimedLight.run_light_step(scene, row, cfg) setup = (
        scene = SIMPLEPLANT.scene;
        cfg = _override_cfg(SIMPLEPLANT.cfg; toricity=true, pixel_size_m=0.40);
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, cfg).rows)
    ) evals = 1

SUITE["Run light"]["step"]["simpleplant scattering links"] =
    @benchmarkable ArchimedLight.run_light_step(scene, row, cfg; scattering_mode=:links) setup = (
        scene = SIMPLEPLANT.scene;
        cfg = _override_cfg(SIMPLEPLANT.cfg; scattering=true, pixel_size_m=0.05);
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, cfg).rows)
    ) evals = 1

SUITE["Run light"]["step"]["sky direct fixture"] =
    @benchmarkable begin
        sky = ArchimedLight.compute_sky(row, cfg)
        turtle = ArchimedLight.build_turtle(cfg, sky)
        ArchimedLight.compute_directional_fluxes(sky, turtle, cfg)
    end setup = (
        cfg = SKY_DIRECT.cfg;
        row = first(SKY_DIRECT.rows)
    ) evals = 1

SUITE["Run light"]["series"]["simpleplant uncached"] =
    @benchmarkable ArchimedLight.run_light_series(scene, meteo, cfg) setup = (
        scene = SIMPLEPLANT.scene;
        meteo = SIMPLEPLANT.meteo;
        cfg = _override_cfg(SIMPLEPLANT.cfg; cache_radiation=false, scattering=false)
    ) evals = 1

SUITE["Run light"]["series"]["simpleplant cached"] =
    @benchmarkable ArchimedLight.run_light_series(scene, meteo, cfg) setup = (
        scene = SIMPLEPLANT.scene;
        meteo = SIMPLEPLANT.meteo;
        cfg = _override_cfg(SIMPLEPLANT.cfg; cache_radiation=true, scattering=false)
    ) evals = 1

SUITE["Synthetic"] = BenchmarkGroup()
SUITE["Synthetic"]["first order stacked plates"] =
    @benchmarkable ArchimedLight.run_light_step(scene, meteo, cfg) setup = (
        scene = SYNTHETIC.scene;
        meteo = SYNTHETIC.meteo;
        cfg = _override_cfg(SYNTHETIC.cfg_ref; sky_sectors=46, all_in_turtle=false, scattering=false, pixel_size_m=0.01)
    ) evals = 1

SUITE["Synthetic"]["stacked plates with scattering"] =
    @benchmarkable ArchimedLight.run_light_step(scene, meteo, cfg; scattering_mode=:raycast) setup = (
        scene = SYNTHETIC.scene;
        meteo = SYNTHETIC.meteo;
        cfg = _override_cfg(SYNTHETIC.cfg_ref; sky_sectors=46, all_in_turtle=false, scattering=true, pixel_size_m=0.01)
    ) evals = 1

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
