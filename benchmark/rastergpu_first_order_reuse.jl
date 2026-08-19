const _BENCH_METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_BENCH_METAL", ""))
const _BENCH_WANTS_METAL = _BENCH_METAL_FLAG in ("1", "true", "yes", "on", "required", "force")

if _BENCH_WANTS_METAL && Base.find_package("Metal") !== nothing
    Base.require(Main, :Metal)
end

using ArchimedLight
using KernelAbstractions
using Printf
using Statistics

function _env_int(name::AbstractString, default::Int)
    return parse(Int, get(ENV, name, string(default)))
end

function _env_float(name::AbstractString, default::Float64)
    return parse(Float64, get(ENV, name, string(default)))
end

const BENCH_SAMPLES = _env_int("ARCHIMEDLIGHT_RASTERGPU_REUSE_SAMPLES", 5)
const BENCH_SECTORS = _env_int("ARCHIMEDLIGHT_BENCH_SECTORS", 6)
const BENCH_PIXEL_SIZE = _env_float("ARCHIMEDLIGHT_BENCH_PIXEL_SIZE", 0.05)
const BENCH_MAX_HITS = _env_int("ARCHIMEDLIGHT_BENCH_MAX_HITS", 32)
const BENCH_TILE_SIZE = _env_int("ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_SIZE", 1)
const BENCH_TILE_FACE_CAPACITY = _env_int("ARCHIMEDLIGHT_BENCH_RASTERGPU_TILE_FACE_CAPACITY", 32)
const BENCH_TOP_HIT_TILE_SIZE = _env_int("ARCHIMEDLIGHT_BENCH_RASTERGPU_TOP_HIT_TILE_SIZE", 1)
const BENCH_TOP_HIT_TILE_FACE_CAPACITY =
    _env_int("ARCHIMEDLIGHT_BENCH_RASTERGPU_TOP_HIT_TILE_FACE_CAPACITY", 32)
const BENCH_OUTPUT = get(ENV, "ARCHIMEDLIGHT_RASTERGPU_REUSE_OUTPUT", "")
const BENCH_WORKLOADS = [
    strip(lowercase(value))
    for value in split(get(ENV, "ARCHIMEDLIGHT_RASTERGPU_REUSE_WORKLOADS", "simple"), ',')
    if !isempty(strip(value))
]
const BENCH_COFFEE_PAVING = _env_int("ARCHIMEDLIGHT_BENCH_COFFEE_PAVING", 25)

function _optional_metal_backend()
    _BENCH_WANTS_METAL || return nothing
    if !Sys.isapple() || Sys.ARCH != :aarch64
        _BENCH_METAL_FLAG in ("required", "force") && error("Metal benchmark requires Apple Silicon macOS.")
        @warn "Metal benchmark requires Apple Silicon macOS; skipping Metal backends."
        return nothing
    end
    if Base.find_package("Metal") === nothing
        _BENCH_METAL_FLAG in ("required", "force") && error("Metal.jl is not available in this environment.")
        @warn "Metal.jl is not available; skipping Metal backends."
        return nothing
    end
    metal_mod = Base.require(Main, :Metal)
    array_type = getproperty(metal_mod, :MtlArray)
    dev_array = Base.invokelatest(array_type, zeros(Float32, 1))
    return Base.invokelatest(KernelAbstractions.get_backend, dev_array)
end

function _options_for(options)
    return ArchimedLight.LightOptions(
        options;
        scattering=false,
        cache_pixel_table=false,
        turtle_sectors=BENCH_SECTORS,
        pixel_size=BENCH_PIXEL_SIZE,
    )
end

function _workload_from_config(name, config; kwargs...)
    options, scene, meteo, models = ArchimedLight.read_config(config; kwargs...)
    options = _options_for(options)
    row = first(ArchimedLight.prepare_meteo(meteo, options).rows)
    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    return (; name, scene, models, options, turtle, fluxes)
end

function _load_workload(name::AbstractString)
    name == "simple" && return _load_simple_workload()
    name == "coffee" && return _load_coffee_workload()
    error("Unsupported ARCHIMEDLIGHT_RASTERGPU_REUSE_WORKLOADS value: $name (supported: simple, coffee)")
end

function _load_simple_workload()
    config = joinpath(@__DIR__, "..", "test", "fast_fixtures", "simpleplant_16_notoric", "input", "config.yml")
    return _workload_from_config("simple", config)
end

function _load_coffee_workload()
    config = joinpath(@__DIR__, "..", "example_2", "config.yml")
    return _workload_from_config("coffee", config; plot_paving_override=BENCH_COFFEE_PAVING)
end

function _time_samples(f, samples::Int)
    times = Float64[]
    last = nothing
    last = f()
    for _ in 1:samples
        elapsed = @elapsed last = f()
        push!(times, elapsed)
    end
    return (; seconds=median(times), min_seconds=minimum(times), max_seconds=maximum(times), result=last)
end

function _scene_metrics(prepared)
    geometry = prepared.geometry
    return (
        nodes=length(geometry.node_ids),
        faces=length(geometry.faces),
        pixels=geometry.plotbox.nx * geometry.plotbox.ny,
    )
end

function _row(workload, name, timing, reference, metrics)
    dense = timing.result.dense
    return (
        workload=workload,
        backend=name,
        nodes=metrics.nodes,
        faces=metrics.faces,
        pixels=metrics.pixels,
        seconds=timing.seconds,
        min_seconds=timing.min_seconds,
        max_seconds=timing.max_seconds,
        hits=sum(dense.hits_per_node),
        hits_equal=reference === nothing ? true : dense.hits_per_node == reference.hits_per_node,
        par_delta=reference === nothing ? 0.0 : maximum(abs.(dense.incident_power.par .- reference.incident_power.par)),
        nir_delta=reference === nothing ? 0.0 : maximum(abs.(dense.incident_power.nir .- reference.incident_power.nir)),
    )
end

function _write_rows(rows)
    isempty(BENCH_OUTPUT) && return rows
    cols = propertynames(first(rows))
    open(BENCH_OUTPUT, "w") do io
        println(io, join(string.(cols), '\t'))
        for row in rows
            println(io, join((string(getproperty(row, col)) for col in cols), '\t'))
        end
    end
    @info "Wrote RasterGPU first-order reuse benchmark" path=BENCH_OUTPUT rows=length(rows)
    return rows
end

function main()
    metal_backend = _optional_metal_backend()
    rows = NamedTuple[]

    for workload_name in BENCH_WORKLOADS
        fixture = _load_workload(workload_name)
        prepared = ArchimedLight._prepare_interception_data(
            fixture.scene,
            fixture.models,
            fixture.options;
            include_budget_maps=true,
        )
        metrics = _scene_metrics(prepared)
        cpu_timing = _time_samples(BENCH_SAMPLES) do
            ArchimedLight._compute_first_order(prepared, fixture.turtle, fixture.fluxes, fixture.options)
        end
        reference = cpu_timing.result.dense
        push!(rows, _row(fixture.name, "raster_cpu_prepared", cpu_timing, nothing, metrics))

        if metal_backend !== nothing
            raster_backend = ArchimedLight.RasterGPUBackend(
                backend=metal_backend,
                max_hits_per_pixel=BENCH_MAX_HITS,
                tile_size=BENCH_TILE_SIZE,
                tile_face_capacity=BENCH_TILE_FACE_CAPACITY,
                top_hit_tile_size=BENCH_TOP_HIT_TILE_SIZE,
                top_hit_tile_face_capacity=BENCH_TOP_HIT_TILE_FACE_CAPACITY,
            )
            raster_data = ArchimedLight._rastergpu_scene_data(prepared, raster_backend.config)
            raster_timing = _time_samples(BENCH_SAMPLES) do
                ArchimedLight.compute_first_order(raster_data, fixture.turtle, fixture.fluxes, fixture.options)
            end
            push!(rows, _row(fixture.name, "rastergpu_metal_prepared", raster_timing, reference, metrics))
        end
    end

    println("RasterGPU first-order prepared-reuse benchmark")
    println("workloads: ", join(BENCH_WORKLOADS, ","))
    println("samples: ", BENCH_SAMPLES)
    println("sectors: ", BENCH_SECTORS)
    println("pixel_size: ", BENCH_PIXEL_SIZE)
    println("tile_size: ", BENCH_TILE_SIZE)
    println("tile_face_capacity: ", BENCH_TILE_FACE_CAPACITY)
    println("top_hit_tile_size: ", BENCH_TOP_HIT_TILE_SIZE)
    println("top_hit_tile_face_capacity: ", BENCH_TOP_HIT_TILE_FACE_CAPACITY)
    println()
    @printf("%-10s %-26s %8s %8s %8s %10s %10s %10s %8s %10s %10s\n",
        "workload", "backend", "nodes", "faces", "pixels", "median", "min", "max", "hits", "PAR d", "NIR d")
    for row in rows
        @printf(
            "%-10s %-26s %8d %8d %8d %10.6f %10.6f %10.6f %8d %10.3e %10.3e\n",
            row.workload,
            row.backend,
            row.nodes,
            row.faces,
            row.pixels,
            row.seconds,
            row.min_seconds,
            row.max_seconds,
            row.hits,
            row.par_delta,
            row.nir_delta,
        )
    end
    return _write_rows(rows)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
