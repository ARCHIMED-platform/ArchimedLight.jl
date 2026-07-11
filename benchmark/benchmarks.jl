using AirspeedVelocity
using ArchimedLight
using BenchmarkTools
using CairoMakie
using Dates
using GeometryBasics
using PlantGeom

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
    rows = collect(ArchimedLight.prepare_meteo(meteo, options))
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
    xs = Float64[]
    ys = Float64[]
    normalized = map(specs) do spec
        p1 = ntuple(j -> Float32(spec.p1[j]), 3)
        p2 = ntuple(j -> Float32(spec.p2[j]), 3)
        p3 = ntuple(j -> Float32(spec.p3[j]), 3)
        p4 = ntuple(j -> Float32(spec.p4[j]), 3)
        append!(xs, (p1[1], p2[1], p3[1], p4[1]))
        append!(ys, (p1[2], p2[2], p3[2], p4[2]))
        mesh = GeometryBasics.Mesh(
            GeometryBasics.Point3f[
                GeometryBasics.Point3f(p1...),
                GeometryBasics.Point3f(p2...),
                GeometryBasics.Point3f(p3...),
                GeometryBasics.Point3f(p4...),
            ],
            GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)],
        )
        (
            mesh=mesh,
            group=String(get(spec, :group, "plate")),
            type=String(get(spec, :type, "plate")),
            object_id=Int(get(spec, :object_id, get(spec, :source_topology_id, 1))),
            source_topology_id=Int(get(spec, :source_topology_id, 1)),
        )
    end

    bounds = (minimum(xs), minimum(ys), maximum(xs), maximum(ys))
    return PlantGeom.make_scene(domain=bounds, source_path="benchmark_synthetic_scene") do builder
        for spec in normalized
            PlantGeom.add_object!(
                builder,
                spec.mesh;
                group=spec.group,
                type=spec.type,
                id=spec.object_id,
                source_topology_id=spec.source_topology_id,
            )
        end
    end
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
        RI_PAR_f=ri_par_f,
        RI_NIR_f=ri_nir_f,
        direct_fraction=direct_fraction,
        sun_azimut=sun_azimut,
        sun_elevation=sun_elevation,
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
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, options))
    ) evals = 1

SUITE["Run light"]["step"]["simpleplant scattering raycast"] =
    @benchmarkable ArchimedLight.run_light_step(scene, models, row, options; scattering_mode=:raycast) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        options = _override_options(SIMPLEPLANT.options; scattering=true, pixel_size_m=0.05);
        row = first(ArchimedLight.prepare_meteo(SIMPLEPLANT.meteo, options))
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

SUITE["Run light"]["series"]["simpleplant prepared cache"] =
    @benchmarkable ArchimedLight.run_light_series(cache, meteo) setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        meteo = SIMPLEPLANT.meteo;
        options = _override_options(SIMPLEPLANT.options; cache_radiation=true, scattering=false);
        cache = ArchimedLight.prepare_light_cache(scene, models, options)
    ) evals = 1

SUITE["Run light"]["series"]["simpleplant prepared cache manual loop"] =
    @benchmarkable begin
        for row in meteo
            ArchimedLight.run_light_step(cache, row)
        end
    end setup = (
        scene = SIMPLEPLANT.scene;
        models = SIMPLEPLANT.models;
        meteo = SIMPLEPLANT.meteo;
        options = _override_options(SIMPLEPLANT.options; cache_radiation=true, scattering=false);
        cache = ArchimedLight.prepare_light_cache(scene, models, options)
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
    @benchmarkable simulate_home_figure(PKG_ROOT) evals = 1

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
